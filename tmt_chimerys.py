#!/usr/bin/env python3

# DIA TMT QUANTIFICATION CHIMERYS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import tomllib
import argparse
import warnings
import pandas as pd
from tqdm import tqdm
from pyteomics import mzml

from typing import Dict
from typing import List
from typing import Tuple
from typing import Any


__version = "0.0.9"
__date = "2025-07-15"

STRATEGY = 1
PROTON = 1.007276466812
TMT = {
    "TMTpro-126": 126.127726,
    "TMTpro-127N": 127.124761,
    "TMTpro-127C": 127.131081,
    "TMTpro-128N": 128.128116,
    "TMTpro-128C": 128.134436,
    "TMTpro-129N": 129.131471,
    "TMTpro-129C": 129.13779,
    "TMTpro-130N": 130.134825,
    "TMTpro-130C": 130.141145,
    "TMTpro-131N": 131.13818,
    "TMTpro-131C": 131.1445,
    "TMTpro-132N": 132.141535,
    "TMTpro-132C": 132.147855,
    "TMTpro-133N": 133.14489,
    "TMTpro-133C": 133.15121,
    "TMTpro-134N": 134.148245,
    "TMTpro-134C": 134.154565,
    "TMTpro-135N": 135.151600,
}
TMT_TOLERANCE = 0.0025


def __get_uncharged_mass_from_exp_mass(mz: float, charge: int) -> float:
    return mz * charge - PROTON * charge


def __read_settings(toml: str) -> Dict[str, Any]:
    parsed_toml = None
    with open(toml, "rb") as f:
        parsed_toml = tomllib.load(f)
        f.close()
    if parsed_toml is None:
        raise RuntimeError("Could not read config file. Is it in valid TOML format?")
    return {
        "window_size": parsed_toml["METHOD"]["window_size"],
        "window_start": parsed_toml["METHOD"]["window_start"],
        "window_end": parsed_toml["METHOD"]["window_end"],
        "mz_tolerance": parsed_toml["MATCHING"]["mz_tolerance"],
        "rt_window": parsed_toml["MATCHING"]["ms1_rt_window"],
        "threshold": parsed_toml["FILTERING"]["total_intensity_threshold"],
        "noise": parsed_toml["FILTERING"]["noise_threshold"],
    }


def __get_tmt_intensities(spectrum: Dict[str, Any]) -> Dict[str, float]:
    # TODO [abundance instead of intensities]
    tmt_quants = {key: 0.0 for key in TMT.keys()}
    for reporter_ion_name, reporter_ion_mass in TMT.items():
        for i, mz in enumerate(spectrum["mz_array"]):
            if __within_tolerance(mz, reporter_ion_mass, TMT_TOLERANCE):
                tmt_quants[reporter_ion_name] += spectrum["intensity_array"][i]
                break
    return tmt_quants


def __get_key(mass: float) -> int:
    return int(round(mass * 10000))


def __parse_scan_nr_from_id(id: str) -> int:
    return int(id.split("scan=")[1].split()[0])


def __read_spectra_by_scannumber(
    filename: str,
) -> Dict[str, Dict[int, Dict[str, Any]]]:
    spectra_ms1 = dict()
    spectra_ms2 = dict()
    total = 0
    total_ms1 = 0
    total_ms2 = 0
    with mzml.read(filename) as reader:
        for spectrum in reader:
            scan = __parse_scan_nr_from_id(spectrum["id"])
            ms_level = int(spectrum["ms level"])
            # check if all fields for retrieving retention time are available
            if (
                "scanList" not in spectrum
                or "scan" not in spectrum["scanList"]
                or len(spectrum["scanList"]["scan"]) != 1
            ):
                raise RuntimeError(f"Can't get retention time for spectrum: {spectrum}")
            # get retention time
            rt_in_min = float(spectrum["scanList"]["scan"][0]["scan start time"])
            # retention time seems to be in minutes in mzML files
            rt_in_sec = rt_in_min * 60.0
            # for MS2 spectra we extract all precursors
            # though realistically there should only be one
            if ms_level == 2:
                if "precursorList" not in spectrum:
                    raise RuntimeError(
                        f"[precursorList] No precursor for MS2 spectrum found: {spectrum}"
                    )
                if (
                    "precursor" not in spectrum["precursorList"]
                    or len(spectrum["precursorList"]["precursor"]) != 1
                ):
                    raise RuntimeError(
                        f"[precursor] No precursor for MS2 spectrum found: {spectrum}"
                    )
                for precursor in spectrum["precursorList"]["precursor"]:
                    if "selectedIonList" not in precursor:
                        raise RuntimeError(
                            f"[selectedIonList] No precursor for MS2 spectrum found: {spectrum}"
                        )
                    if (
                        "selectedIon" not in precursor["selectedIonList"]
                        or len(precursor["selectedIonList"]["selectedIon"]) != 1
                    ):
                        raise RuntimeError(
                            f"[selectedIon] No precursor for MS2 spectrum found: {spectrum}"
                        )
                    for ion in precursor["selectedIonList"]["selectedIon"]:
                        s = dict()
                        s["scan_nr"] = scan
                        s["precursor"] = float(ion["selected ion m/z"])
                        s["rt"] = rt_in_sec
                        s["mz_array"] = spectrum["m/z array"]
                        s["intensity_array"] = spectrum["intensity array"]
                        if scan in spectra_ms2:
                            raise RuntimeError()
                        else:
                            spectra_ms2[scan] = s
                            total += 1
                            total_ms2 += 1

            # for MS1 spectra we save them based on retention time
            elif ms_level == 1:
                # the primary key is the retention time in seconds * 10 000 (rounded)
                primary_key = __get_key(rt_in_sec)
                s = dict()
                s["scan_nr"] = scan
                s["precursor"] = None
                s["rt"] = rt_in_sec
                s["mz_array"] = spectrum["m/z array"]
                s["intensity_array"] = spectrum["intensity array"]
                # if there are multiple MS1 spectra with the same retention time
                # an error is raised
                if primary_key in spectra_ms1:
                    raise KeyError(
                        f"MS1 spectrum for retention time {s['rt']} already exists!"
                    )
                else:
                    spectra_ms1[primary_key] = s
                    total += 1
                    total_ms1 += 1
            else:
                raise ValueError(
                    f"Found spectrum with MS level {ms_level}. Not supported!"
                )
    print(f"Total number of parsed spectra: {total}")
    print(f"Number of MS1 spectra: {total_ms1}")
    print(f"Number of MS2 spectra: {total_ms2}")
    return {"ms1": spectra_ms1, "ms2": spectra_ms2}


# checks if two values (value and ref) are equal with a certain tolerance (tol)
def __within_tolerance(value: float, ref: float, tol: float) -> bool:
    return (value > ref - tol) and (value < ref + tol)


# checks if the specified precursor m/z can be found in the given MS1 spectrum
# considering tolerance
def __check_mz_in_ms1(
    precursor_mz: float, spectrum: Dict[str, Any], mz_tol: float
) -> bool:
    for mz in spectrum["mz_array"]:
        if __within_tolerance(mz, precursor_mz, mz_tol):
            return True
    return False


# calculates the purity for a given precursor
def __calculate_precursor_intensity_ms1(
    precursor_mz: float,
    spectrum: Dict[str, Any],  # MS1 spectrum
    mz_tol: float,
    noise_threshold: float,
    windows: List[Tuple[float, float]],
) -> float:
    # parameter windows should define all DIA windows including their start and
    # end points (in m/z)
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
    matching_window = None
    for window in windows:
        if precursor > window[0] and precursor <= window[1]:
            matching_window = window
            break
    # if no matching window is found an error is raised
    if matching_window is None:
        raise RuntimeError(
            f"Could not find matching window for precursor m/z {precursor_mz}!"
        )
    # get highest intensity peak in window
    most_intense_peak = 0.0
    for i in range(len(spectrum["mz_array"])):
        if (
            spectrum["mz_array"][i] > matching_window[0]
            and spectrum["mz_array"][i] <= matching_window[1]
        ):
            if spectrum["intensity_array"][i] > most_intense_peak:
                most_intense_peak = spectrum["intensity_array"][i]
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
    noise_threshold: float,
    spectra: Dict[str, Any],
    windows: List[Tuple[float, float]],
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
        precursor_mz, ms1, mz_tol, noise_threshold, windows
    )
    return {"spectrum": spectrum, "purity": purity}


# given window start and end points (in m/z) and a window size, calculate all
# m/z windows
def __get_windows(
    window_start: float, window_end: float, window_size: float
) -> List[Tuple[float, float]]:
    windows = list()
    current_window_start = window_start
    while current_window_start < window_end:
        current_window_end = current_window_start + window_size
        if current_window_end > window_end:
            windows.append((current_window_start, window_end))
            break
        windows.append((current_window_start, current_window_end))
        current_window_start += window_size
    return windows


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
        # get all m/z windows
        windows = __get_windows(
            float(settings["window_start"]),
            float(settings["window_end"]),
            float(settings["window_size"]),
        )
        # get corresponding MS2 spectrum for identification
        spectrum_purity = __get_ms2_spectrum_by_scannumber(
            scan_nr,
            prec_mz,
            mz_tol,
            rt_window,
            noise_threshold,
            spectra,
            windows,
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
        else:
            if quantify:
                for key in channels.keys():
                    channels[key].append(None)
            purities.append(None)
            nr_of_missing_ms1 += 1
    # update Chimerys result
    df["Co-Isolation Purity"] = purities
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
        prog="tmt_chimerys.py",
        description="Calculates co-isolation purity for Chimerys DIA TMT PSMs and optionally quantifies them.",
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
