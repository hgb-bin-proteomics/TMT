#!/usr/bin/env python3

# UNNAMED TMT PROJECT
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


__version = "0.0.5"
__date = "2025-07-09"

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
        "precursor_mass": parsed_toml["SPECTRONAUT"]["precursor_mass"],
        "precursor_mz": parsed_toml["SPECTRONAUT"]["precursor_mz"],
        "precursor_charge": parsed_toml["SPECTRONAUT"]["precursor_mz"],
        "retention_time": parsed_toml["SPECTRONAUT"]["retention_time"],
        "retention_time_in_sec": parsed_toml["SPECTRONAUT"]["retention_time_in_sec"],
        "mz_tolerance": parsed_toml["MATCHING"]["mz_tolerance"],
        "rt_tolerance": parsed_toml["MATCHING"]["rt_tolerance"],
        "rt_window": parsed_toml["MATCHING"]["ms1_rt_window"],
        "threshold": parsed_toml["FILTERING"]["total_intensity_threshold"],
        "noise": parsed_toml["FILTERING"]["noise_threshold"],
    }


def __get_tmt_intensities(spectrum: Dict[str, Any]) -> Dict[str, float]:
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
                if "precursorList" in spectrum:
                    for precursor in spectrum["precursorList"]["precursor"]:
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
                else:
                    # an error is raised if no precursor is found
                    raise KeyError(f"No precursor for MS2 spectrum found: {spectrum}")
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


# check the intensity of the MS2 precursor in the MS1 spectrum if it passes the
# specified threshold -> where the threshold is defined as precursor intensity
# divided by total intensity in the window
def __check_precursor_intensity_ms1(
    precursor_mz: float,
    spectrum: Dict[str, Any],
    mz_tol: float,
    filter_threshold: float,
    noise_threshold: float,
    windows: List[Tuple[float, float]],
) -> bool:
    # parameter windows should define all DIA windows including their start and
    # end points (in m/z)
    precursor_index = None
    precursor_intensity = None
    # first look for the precursor in the MS1 spectrum
    for i in range(len(spectrum["mz_array"])):
        if __within_tolerance(spectrum["mz_array"][i], precursor_mz, mz_tol):
            # if no precursor yet found, set precursor
            if precursor_index is None:
                precursor_index = i
                precursor_intensity = spectrum["intensity_array"][i]
            # else we need to deal with the fact that there is multiple potential
            # precursors
            else:
                # todo clarify behaviour
                if STRATEGY == 1:
                    # use precursor with highest intensity
                    if spectrum["intensity_array"][i] > precursor_intensity:
                        precursor_index = i
                        precursor_intensity = spectrum["intensity_array"][i]
                elif STRATEGY == 2:
                    # do not use identification
                    return False
                elif STRATEGY == 3:
                    # use closest precursor
                    raise NotImplementedError()
                else:
                    raise RuntimeError(
                        f"Found ambiguous precursors in MS1 spectrum for precursor m/z {precursor_mz} using MS1 spectrum at retention time {spectrum['rt']}."
                    )
    # if no precursor is found an error is raised
    if precursor_index is None or precursor_intensity is None:
        raise RuntimeError("Could not find precursor in MS1 spectrum.")
    # find the corresponding m/z window that the precursor is in
    matching_window = None
    for window in windows:
        if precursor_mz > window[0] and precursor_mz < window[1]:
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
            and spectrum["mz_array"][i] < matching_window[1]
        ):
            if spectrum["intensity_array"][i] > most_intense_peak:
                most_intense_peak = spectrum["intensity_array"][i]
    # if precursor is noisy, return False
    if precursor_intensity / most_intense_peak < noise_threshold:
        return False
    # calculate total intensity in window
    total_intensity_in_window = 0.0
    for i in range(len(spectrum["mz_array"])):
        if (
            spectrum["mz_array"][i] > matching_window[0]
            and spectrum["mz_array"][i] < matching_window[1]
        ):
            if spectrum["intensity_array"][i] / most_intense_peak > noise_threshold:
                total_intensity_in_window += spectrum["intensity_array"][i]
    # return if intensity ratio passes threshold
    return precursor_intensity / total_intensity_in_window > filter_threshold


# get corresponding MS2 spectrum for a Spectronaut given precursor m/z and
# retention time
# May return None if no MS2 spectrum is found or if it does not pass the
# intensity filter
def __get_ms2_spectrum(
    precursor_mz: float,
    mz_tol: float,
    retention_time: float,
    rt_tol: float,
    retention_time_ms1_window: float,
    window_size_unidirectional: float,
    filter_threshold: float,
    noise_threshold: float,
    spectra: Dict[str, Any],
    windows: List[Tuple[float, float]],
    verbose: int,
) -> Dict[str, Any] | None:
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
        if verbose == 2:
            raise RuntimeError(
                f"Could not find a suitable precursor for precursor m/z {precursor_mz} and retention time {retention_time}."
            )
        elif verbose == 1:
            warnings.warn(
                RuntimeWarning(
                    f"Could not find a suitable precursor for precursor m/z {precursor_mz} and retention time {retention_time}."
                )
            )
            return None
        else:
            return None
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
        if verbose == 2:
            raise RuntimeError(
                f"Could not find a suitable retention time for precursor m/z {precursor_mz} and retention time {retention_time}."
            )
        elif verbose == 1:
            warnings.warn(
                RuntimeWarning(
                    f"Could not find a suitable retention time for precursor m/z {precursor_mz} and retention time {retention_time}."
                )
            )
            return None
        else:
            return None
    # look for closest MS1 spectrum that has precursor (with m/z tolerance mz_tol)
    # within ms1 rt window
    retention_time_ms1_window_range = __get_key(retention_time_ms1_window)
    ms1 = None
    for i in range(retention_time_ms1_window_range):
        if secondary_key_base - i in spectra["ms1"]:
            if __check_mz_in_ms1(
                precursor_mz, spectra["ms1"][secondary_key_base - i], mz_tol
            ):
                ms1 = spectra["ms1"][secondary_key_base - i]
                break
        if secondary_key_base + i in spectra["ms1"]:
            if __check_mz_in_ms1(
                precursor_mz, spectra["ms1"][secondary_key_base + i], mz_tol
            ):
                ms1 = spectra["ms1"][secondary_key_base + i]
                break
    # if no MS1 spectrum is found -> no error, but None returned
    if ms1 is None:
        warnings.warn(
            RuntimeWarning(
                f"Could not find a suitable MS1 spectrum for precursor m/z {precursor_mz} and retention time {retention_time}."
            )
        )
        return None
    # intensity filter
    if __check_precursor_intensity_ms1(
        precursor_mz, ms1, mz_tol, filter_threshold, noise_threshold, windows
    ):
        return spectrum
    return None


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
    nr_of_filtered_spectra = 0
    for i, row in tqdm(
        df.iterrows(), total=df.shape[0], desc="Annotating Spectronaut result..."
    ):
        # get Spectronaut precursor m/z from identification
        prec_mz = float(row[settings["precursor_mz"]])
        # settings defined m/z tolerance
        mz_tol = float(settings["mz_tolerance"])
        # get Spectronaut retention time from identification
        rt = float(row[settings["retention_time"]])
        # settings  defined m/z tolerance in seconds
        rt_tol = float(settings["rt_tolerance"])
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
        # convert retention time to seconds
        if not settings["retention_time_in_sec"]:
            rt = rt * 60.0
        # get corresponding MS2 spectrum for identification
        spectrum = __get_ms2_spectrum(
            prec_mz,
            mz_tol,
            rt,
            rt_tol,
            rt_window,
            window_size_unidirectional,
            filter_threshold,
            noise_threshold,
            spectra,
            windows,
            verbose,
        )
        # if a spectrum that passes the intensity filter is found -> quantify
        if spectrum is not None:
            tmt_quants = __get_tmt_intensities(spectrum)
            for key in channels.keys():
                channels[key].append(tmt_quants[key])
        else:
            for key in channels.keys():
                channels[key].append(None)
            nr_of_filtered_spectra += 1
    # update Spectronaut result
    for key in channels.keys():
        df[key] = channels[key]
    print(
        f"Total number of identifications with MS1 spectra below threshold: {nr_of_filtered_spectra}"
    )
    return df


def main(argv=None) -> pd.DataFrame:
    parser = argparse.ArgumentParser(
        prog="TMT",
        description="TMT",
        epilog="Bottom Text",
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
        "-v",
        "--verbose",
        dest="verbose",
        default=2,
        help="Verbose level.",
        type=int,
    )
    parser.add_argument(
        "-w",
        "--window",
        dest="window",
        default=None,
        help="Window size, overrides config file!",
        type=float,
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
    df.to_csv(args.spectronaut + "_tmt_quant.csv", index=False)
    return df


if __name__ == "__main__":
    _ = main()
