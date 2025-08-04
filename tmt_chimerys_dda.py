#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "pandas",
#   "tqdm",
#   "pyteomics[XML]",
# ]
# ///

# DDA TMT QUANTIFICATION CHIMERYS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import argparse
import warnings
import pandas as pd
from tqdm import tqdm

from typing import Dict
from typing import Any

from tmt_chimerys import TMT
from tmt_chimerys import __get_bool_from_value
from tmt_chimerys import __read_settings
from tmt_chimerys import __get_tmt_intensities
from tmt_chimerys import __read_spectra_by_scannumber
from tmt_chimerys import __get_key
from tmt_chimerys import __within_tolerance
from tmt_chimerys import __check_mz_in_ms1

__version = "0.0.10"
__date = "2025-07-22"

ISOTOPE = 1.00335
STRATEGY = 1


# calculates the purity for a given precursor
def __calculate_precursor_intensity_ms1(
    precursor_mz: float,
    spectrum: Dict[str, Any],  # MS1 spectrum
    mz_tol: float,
    do_deisotope: bool,
    isotope_tol: float,
    max_charge: int,
    noise_threshold: float,
    window_size: float,
) -> float:
    precursor = None
    precursor_index = None
    precursor_intensity = None
    # first look for the precursor in the MS1 spectrum
    for i in range(len(spectrum["mz_array"])):
        if __within_tolerance(spectrum["mz_array"][i], precursor_mz, mz_tol):
            # if no precursor yet found, set precursor
            if precursor is None:
                precursor = spectrum["mz_array"][i]
                precursor_index = i
                precursor_intensity = spectrum["intensity_array"][i]
            # else we need to deal with the fact that there is multiple potential
            # precursors
            else:
                # todo clarify behaviour
                if STRATEGY == 1:
                    # use precursor with highest intensity
                    if spectrum["intensity_array"][i] > precursor_intensity:
                        precursor = spectrum["mz_array"][i]
                        precursor_index = i
                        precursor_intensity = spectrum["intensity_array"][i]
                elif STRATEGY == 2:
                    # do not use identification
                    return -1.0
                elif STRATEGY == 3:
                    # use closest precursor
                    raise NotImplementedError()
                else:
                    raise RuntimeError(
                        f"Found ambiguous precursors in MS1 spectrum for precursor m/z {precursor_mz} using MS1 spectrum at retention time {spectrum['rt']}."
                    )
    # if no precursor is found an error is raised
    if precursor is None or precursor_index is None or precursor_intensity is None:
        raise RuntimeError("Could not find precursor in MS1 spectrum.")
    # find the corresponding m/z window that the precursor is in
    matching_window = (precursor - window_size / 2.0, precursor + window_size / 2.0)
    # get highest intensity peak in window
    most_intense_peak = 0.0
    peaks_in_window_mz = list()
    peaks_in_window_in = list()
    for i in range(len(spectrum["mz_array"])):
        if (
            spectrum["mz_array"][i] > matching_window[0]
            and spectrum["mz_array"][i] <= matching_window[1]
        ):
            if spectrum["intensity_array"][i] > most_intense_peak:
                most_intense_peak = spectrum["intensity_array"][i]
            peaks_in_window_mz.append(spectrum["mz_array"][i])
            peaks_in_window_in.append(spectrum["intensity_array"][i])
    if do_deisotope:
        # look for precursor isotope peaks
        charge = None
        ## determine charge
        for i in range(1, 11):  # +- 10 peaks are considered for charge determination
            if precursor_index + i < len(spectrum["mz_array"]):
                curr_peak_mz = spectrum["mz_array"][precursor_index + i]
                for c in range(1, max_charge + 1):
                    if __within_tolerance(
                        curr_peak_mz, precursor + ISOTOPE / c, isotope_tol
                    ):
                        charge = c
                        break
            if charge is None:
                if precursor_index - i >= 0:
                    curr_peak_mz = spectrum["mz_array"][precursor_index - i]
                    for c in range(1, max_charge + 1):
                        if __within_tolerance(
                            curr_peak_mz, precursor - ISOTOPE / c, isotope_tol
                        ):
                            charge = c
                            break
            if charge is not None:
                break
        # get precursor isotope peaks
        if charge is not None:
            for i in range(1, 6):  # +- 5 isotopes are considered for every precursor
                precursor_isotope_plus_i = precursor + i * (ISOTOPE / charge)
                if precursor_isotope_plus_i <= matching_window[1]:
                    for j in range(len(peaks_in_window_mz)):
                        if peaks_in_window_in[j] / most_intense_peak >= noise_threshold:
                            if __within_tolerance(
                                peaks_in_window_mz[j],
                                precursor_isotope_plus_i,
                                isotope_tol,
                            ):
                                precursor_intensity += peaks_in_window_in[j]
                                break
                precursor_isotope_minus_i = precursor - i * (ISOTOPE / charge)
                if precursor_isotope_minus_i > matching_window[0]:
                    for j in range(len(peaks_in_window_mz)):
                        if peaks_in_window_in[j] / most_intense_peak >= noise_threshold:
                            if __within_tolerance(
                                peaks_in_window_mz[j],
                                precursor_isotope_minus_i,
                                isotope_tol,
                            ):
                                precursor_intensity += peaks_in_window_in[j]
                                break
    # if precursor is noisy, return 0.0
    if precursor_intensity / most_intense_peak < noise_threshold:
        return 0.0
    # calculate total intensity in window
    total_intensity_in_window = 0.0
    for i in range(len(spectrum["mz_array"])):
        if (
            spectrum["mz_array"][i] > matching_window[0]
            and spectrum["mz_array"][i] <= matching_window[1]
        ):
            if spectrum["intensity_array"][i] / most_intense_peak > noise_threshold:
                total_intensity_in_window += spectrum["intensity_array"][i]
    # return intensity ratio (purity)
    return precursor_intensity / total_intensity_in_window


def __get_ms2_spectrum_by_scannumber(
    scan_nr: int,
    precursor_mz: float,
    mz_tol: float,
    retention_time_ms1_window: float,
    do_deisotope: bool,
    isotope_tol: float,
    max_charge: int,
    noise_threshold: float,
    spectra: Dict[str, Any],
    window_size: float,
) -> Dict[str, Any]:
    spectrum = spectra["ms2"][scan_nr]
    retention_time_key = __get_key(spectrum["rt"])
    # look for closest MS1 spectrum that has precursor (with m/z tolerance mz_tol)
    # within ms1 rt window
    retention_time_ms1_window_range = __get_key(retention_time_ms1_window)
    ms1 = None
    for i in range(retention_time_ms1_window_range):
        if retention_time_key - i in spectra["ms1"]:
            if __check_mz_in_ms1(
                precursor_mz, spectra["ms1"][retention_time_key - i], mz_tol
            ):
                ms1 = spectra["ms1"][retention_time_key - i]
                break
        if retention_time_key + i in spectra["ms1"]:
            if __check_mz_in_ms1(
                precursor_mz, spectra["ms1"][retention_time_key + i], mz_tol
            ):
                ms1 = spectra["ms1"][retention_time_key + i]
                break
    # if no MS1 spectrum is found -> no error, but None returned
    if ms1 is None:
        warnings.warn(
            RuntimeWarning(
                f"Could not find a suitable MS1 spectrum for precursor m/z {precursor_mz} and retention time {spectrum['rt']}. MS2 scan: {scan_nr}"
            )
        )
        return {"spectrum": None, "purity": None}
    # intensity filter
    purity = __calculate_precursor_intensity_ms1(
        precursor_mz,
        ms1,
        mz_tol,
        do_deisotope,
        isotope_tol,
        max_charge,
        noise_threshold,
        window_size,
    )
    return {"spectrum": spectrum, "purity": purity}


# annotates the Chimerys result with purity and TMT quantities
# currently based on the PSM table
def __annotate_chimerys_result(
    filename: str,
    spectra: Dict[str, Any],
    settings: Dict[str, Any],
    quantify: bool = False,
) -> pd.DataFrame:
    if quantify:
        print("Quantification enabled!")
    else:
        print("Quantification disabled! Only calculating co-isolation purity!")
    # spectra should be given by __read_spectra_by_scannumber
    # settings should be given by __read_settings
    df = pd.read_csv(filename, sep="\t", low_memory=False)
    channels = {key: [] for key in TMT.keys()}
    purities = list()
    parsed_scannr = list()
    nr_of_missing_ms1 = 0
    nr_of_impure_ids = 0
    for i, row in tqdm(
        df.iterrows(), total=df.shape[0], desc="Annotating Chimerys result..."
    ):
        # get scan number
        scan_nr = int(row["First Scan"])
        # get Chimerys precursor m/z from identification
        prec_mz = float(row["m/z [Da]"])
        # settings defined m/z tolerance
        mz_tol = float(settings["mz_tolerance"])
        # settings definied retention time window for MS1 spectra
        rt_window = float(settings["rt_window"])
        # settings defined precursor intensity ratio threshold
        filter_threshold = float(settings["threshold"])
        # settings defined noise threshold
        noise_threshold = float(settings["noise"])
        # get window size
        window_size = float(settings["window_size"])
        # isotope parameters
        do_deisotope = __get_bool_from_value(settings["deisotope"])
        isotope_tolerance = float(settings["isotope_tolerance"])
        max_charge = int(settings["max_charge"])
        # get corresponding MS2 spectrum for identification
        spectrum_purity = __get_ms2_spectrum_by_scannumber(
            scan_nr,
            prec_mz,
            mz_tol,
            rt_window,
            do_deisotope,
            isotope_tolerance,
            max_charge,
            noise_threshold,
            spectra,
            window_size,
        )
        spectrum = spectrum_purity["spectrum"]
        purity = spectrum_purity["purity"]
        # if a spectrum is found -> purity and quantify
        if spectrum is not None:
            if quantify:
                tmt_quants = __get_tmt_intensities(spectrum)
                for key in channels.keys():
                    channels[key].append(tmt_quants[key])
            purities.append(purity)
            if purity < filter_threshold:
                nr_of_impure_ids += 1
            parsed_scannr.append(spectrum["scan_nr"])
        else:
            if quantify:
                for key in channels.keys():
                    channels[key].append(None)
            purities.append(None)
            nr_of_missing_ms1 += 1
            parsed_scannr.append(None)
    # update Chimerys result
    df["Co-Isolation Purity"] = purities
    df["Parsed MS2 Scan Number"] = parsed_scannr
    if quantify:
        for key in channels.keys():
            df[f"Annotated {key}"] = channels[key]
    print(f"Total number of identifications: {df.shape[0]}")
    print(f"Total number of identifications with impure precursors: {nr_of_impure_ids}")
    print(
        f"Total number of identifications with MS1 spectra not found: {nr_of_missing_ms1}"
    )
    return df


def main(argv=None) -> pd.DataFrame:
    parser = argparse.ArgumentParser(
        prog="tmt_chimerys_dda.py",
        description="Calculates co-isolation purity for Chimerys DDA TMT PSMs and optionally quantifies them.",
        epilog="(c) Research Institute of Molecular Pathology, 2025",
    )
    parser.add_argument(
        "-i",
        "--chimerys",
        dest="chimerys",
        required=True,
        help="Path/name of the Chimerys result file in tab-separated .txt format.",
        type=str,
    )
    parser.add_argument(
        "-s",
        "--spectra",
        dest="spectra",
        required=True,
        help="Path/name of the mass spectra file in mzML format.",
        type=str,
    )
    parser.add_argument(
        "-c",
        "--config",
        dest="config",
        required=True,
        help="Path/name of the config file.",
        type=str,
    )
    parser.add_argument(
        "-w",
        "--window",
        dest="window",
        default=None,
        help="Window size, overrides config file!",
        type=float,
    )
    parser.add_argument(
        "-q",
        "--quantify",
        dest="quantify",
        default=False,
        action="store_true",
        help="Enables separate TMT quantification.",
    )
    parser.add_argument("--version", action="version", version=__version)
    args = parser.parse_args(argv)
    settings = __read_settings(args.config)
    if args.window is not None:
        settings["window_size"] = float(args.window)
    print(settings)
    spectra = __read_spectra_by_scannumber(args.spectra)
    df = __annotate_chimerys_result(args.chimerys, spectra, settings, args.quantify)
    df.to_csv(
        args.chimerys.split(".txt")[0] + "_purity_tmt_quant.txt",
        sep="\t",
        index=False,
    )
    return df


if __name__ == "__main__":
    _ = main()
