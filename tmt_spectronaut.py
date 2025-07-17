#!/usr/bin/env python3

# DIA TMT QUANTIFICATION SPECTRONAUT
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import argparse
import warnings
import pandas as pd
from tqdm import tqdm
from pyteomics import mzml

from typing import Dict
from typing import List
from typing import Tuple
from typing import Any

from tmt_chimerys import TMT
from tmt_chimerys import __read_settings
from tmt_chimerys import __get_tmt_intensities
from tmt_chimerys import __get_key
from tmt_chimerys import __check_mz_in_ms1
from tmt_chimerys import __calculate_precursor_intensity_ms1
from tmt_chimerys import __get_windows

__version = "0.0.9"
__date = "2025-07-17"

RT_TOL = 3.0


# read mass spectra from an mzML file
def __read_spectra(filename: str) -> Dict[str, Any]:
    # reads an mzML file and sorts spectra into an easily accessible data structure
    # that stores them based on precursor m/z and retention time
    # returned is a Dict[str, Dict[int, Any]]
    # at the top level there are two keys: "ms1" whichs maps to a dictionary containing
    # all MS1 spectra
    # and "ms2" which maps to a dictionary containing all MS2 spectra
    # MS1: Dict[int, Spectrum]
    # maps retention time in seconds * 10 000 (rounded) to MS1 spectra
    # MS2: Dict[int, Dict[int, Spectrum]]
    # maps precursor m/z * 10 000 (rounded) to dictionaries that map retention time
    # in seconds * 10 000 (rounded) to MS2 spectra
    spectra_ms1 = dict()
    spectra_ms2 = dict()
    total = 0
    total_ms1 = 0
    total_ms2 = 0
    with mzml.read(filename) as reader:
        for spectrum in reader:
            # get MS level
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
                        # the primary key is the precursor m/z * 10 000 (rounded)
                        primary_key = __get_key(float(ion["selected ion m/z"]))
                        # the secondary key is the retention time in seconds * 10 000 (rounded)
                        secondary_key = __get_key(rt_in_sec)
                        s = dict()
                        s["precursor"] = float(ion["selected ion m/z"])
                        s["rt"] = rt_in_sec
                        s["mz_array"] = spectrum["m/z array"]
                        s["intensity_array"] = spectrum["intensity array"]
                        if primary_key in spectra_ms2:
                            # if there already is another MS2 spectrum with the same precursor m/z
                            # and same retention time an error is raised
                            if secondary_key in spectra_ms2[primary_key]:
                                raise KeyError(
                                    f"MS2 spectrum for precursor {s['precursor']} and retention time {s['rt']} already exists!"
                                )
                            else:
                                spectra_ms2[primary_key][secondary_key] = s
                                total += 1
                                total_ms2 += 1
                        else:
                            spectra_ms2[primary_key] = {secondary_key: s}
                            total += 1
                            total_ms2 += 1
            # for MS1 spectra we save them based on retention time
            elif ms_level == 1:
                # the primary key is the retention time in seconds * 10 000 (rounded)
                primary_key = __get_key(rt_in_sec)
                s = dict()
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


# get corresponding MS2 spectrum for a Spectronaut given precursor m/z and
# retention time
# May return None if no MS2 spectrum or corresponding MS1 spectrum is found
def __get_ms2_spectrum(
    precursor_mz: float,
    mz_tol: float,
    retention_time: float,
    rt_tol: float,
    retention_time_ms1_window: float,
    window_size_unidirectional: float,
    noise_threshold: float,
    spectra: Dict[str, Any],
    windows: List[Tuple[float, float]],
    verbose: int,
) -> Dict[str, Any]:
    # get by most similar precursor m/z within window
    primary_key_base = __get_key(precursor_mz)
    primary_key_window = __get_key(window_size_unidirectional)
    precursor = None
    # take the closest precursor m/z MS2 spectrum within the window_size (half)
    for i in range(primary_key_window):
        if primary_key_base - i in spectra["ms2"]:
            precursor = spectra["ms2"][primary_key_base - i]
            break
        if primary_key_base + i in spectra["ms2"]:
            precursor = spectra["ms2"][primary_key_base + i]
            break
    # if no precursor is found, raise an error
    if precursor is None:
        if verbose > 1:
            raise RuntimeError(
                f"[no precursor] Could not find a suitable precursor for precursor m/z {precursor_mz} and retention time {retention_time}."
            )
        elif verbose == 1:
            warnings.warn(
                RuntimeWarning(
                    f"[no precursor] Could not find a suitable precursor for precursor m/z {precursor_mz} and retention time {retention_time}."
                )
            )
            return {"spectrum": None, "purity": None}
        else:
            return {"spectrum": None, "purity": None}
    # get by most similar retention time using a retention time tolerance
    secondary_key_base = __get_key(retention_time)
    secondary_key_window = __get_key(rt_tol)
    spectrum = None
    # take the MS2 spectrum closest in retention time but maximum rt_tol seconds away
    for i in range(secondary_key_window):
        if secondary_key_base - i in precursor:
            spectrum = precursor[secondary_key_base - i]
            break
        if secondary_key_base + i in precursor:
            spectrum = precursor[secondary_key_base + i]
            break
    # if no MS2 spectrum is found an error is raised
    if spectrum is None:
        if verbose > 1:
            raise RuntimeError(
                f"[no rt] Could not find a suitable retention time for precursor m/z {precursor_mz} and retention time {retention_time}."
            )
        elif verbose == 1:
            warnings.warn(
                RuntimeWarning(
                    f"[no rt] Could not find a suitable retention time for precursor m/z {precursor_mz} and retention time {retention_time}."
                )
            )
            return {"spectrum": None, "purity": None}
        else:
            return {"spectrum": None, "purity": None}
    # look for closest MS1 spectrum that has precursor (with m/z tolerance mz_tol)
    # within ms1 rt window
    retention_time_key = __get_key(spectrum["rt"])
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
                f"Could not find a suitable MS1 spectrum for precursor m/z {precursor_mz} and retention time {spectrum['rt']}."
            )
        )
        return {"spectrum": None, "purity": None}
    # intensity filter
    purity = __calculate_precursor_intensity_ms1(
        precursor_mz, ms1, mz_tol, noise_threshold, windows
    )
    return {"spectrum": spectrum, "purity": purity}


# annotates the Spectronaut result with TMT quantities
# currently based on the factory report
def __annotate_spectronaut_result(
    spectronaut_filename: str,
    spectra: Dict[str, Any],
    settings: Dict[str, Any],
    verbose: int = 2,
) -> pd.DataFrame:
    # spectra should be given by __read_spectra
    # settings should be given by __read_settings
    df = pd.read_csv(spectronaut_filename, low_memory=False)
    channels = {key: [] for key in TMT.keys()}
    purities = list()
    nr_of_missing_ms1 = 0
    nr_of_impure_ids = 0
    for i, row in tqdm(
        df.iterrows(), total=df.shape[0], desc="Annotating Spectronaut result..."
    ):
        # get Spectronaut precursor m/z from identification
        prec_mz = float(row["FG.PrecMz"])
        # settings defined m/z tolerance
        mz_tol = float(settings["mz_tolerance"])
        # get Spectronaut retention time from identification
        rt = float(row["EG.ApexRT"]) * 60.0
        # settings  defined m/z tolerance in seconds
        rt_tol = RT_TOL
        # settings definied retention time window for MS1 spectra
        rt_window = float(settings["rt_window"])
        # settings defined DIA window size
        window_size_unidirectional = float(settings["window_size"]) / 2.0
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
        spectrum_purity = __get_ms2_spectrum(
            prec_mz,
            mz_tol,
            rt,
            rt_tol,
            rt_window,
            window_size_unidirectional,
            noise_threshold,
            spectra,
            windows,
            verbose,
        )
        spectrum = spectrum_purity["spectrum"]
        purity = spectrum_purity["purity"]
        # if a spectrum is found -> purity and quantify
        if spectrum is not None:
            tmt_quants = __get_tmt_intensities(spectrum)
            for key in channels.keys():
                channels[key].append(tmt_quants[key])
            purities.append(purity)
            if purity < filter_threshold:
                nr_of_impure_ids += 1
        else:
            for key in channels.keys():
                channels[key].append(None)
            purities.append(None)
            nr_of_missing_ms1 += 1
    # update Spectronaut result
    df["Co-Isolation Purity"] = purities
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
        prog="tmt_spectronaut.py",
        description="Calculates co-isolation purity for Spectronaut DIA TMT peptide matches and quantifies them.",
        epilog="(c) Research Institute of Molecular Pathology, 2025",
    )
    parser.add_argument(
        "-i",
        "--spectronaut",
        dest="spectronaut",
        required=True,
        help="Path/name of the Spectronaut result file.",
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
        default=None,
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
        "-v",
        "--verbose",
        dest="verbose",
        default=2,
        help="Verbose level.",
        type=int,
    )
    parser.add_argument("--version", action="version", version=__version)
    args = parser.parse_args(argv)
    settings = __read_settings(args.config)
    if args.window is not None:
        settings["window_size"] = float(args.window)
    print(settings)
    spectra = __read_spectra(args.spectra)
    df = __annotate_spectronaut_result(
        args.spectronaut, spectra, settings, int(args.verbose)
    )
    df.to_csv(args.spectronaut.split(".csv")[0] + "_purity_tmt_quant.csv", index=False)
    return df


if __name__ == "__main__":
    _ = main()
