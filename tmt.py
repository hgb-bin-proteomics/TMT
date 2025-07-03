#!/usr/bin/env python3

# UNNAMED TMT PROJECT
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import tomllib
import argparse
import pandas as pd
from tqdm import tqdm
from pyteomics import mzml

from typing import Dict
from typing import Any


__version = "0.0.2"
__date = "2025-07-02"

F = "20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz5_3ng_1 1.mzML"
S = "20250613_125208_TT_multi003_mz5_rep1_Birkl_Factory_Report.csv"
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


def __get_uncharged_mass_from_exp_mass(mz: float, charge: int) -> float:
    return mz * charge - PROTON * charge


def __read_settings(toml: str) -> Dict[str, Any]:
    parsed_toml = None
    with open(toml, "rb") as f:
        parsed_toml = tomllib.load(f)
        f.close()
    if parsed_toml is None:
        raise RuntimeError()
    return {
        "window_size": parsed_toml["METHOD"]["window"],
        "precursor_mass": parsed_toml["SPECTRONAUT"]["precursor_mass"],
        "precursor_mz": parsed_toml["SPECTRONAUT"]["precursor_mz"],
        "precursor_charge": parsed_toml["SPECTRONAUT"]["precursor_mz"],
        "retention_time": parsed_toml["SPECTRONAUT"]["retention_time"],
        "retention_time_in_sec": parsed_toml["SPECTRONAUT"]["retention_time_in_sec"],
        "mz_tolerance": parsed_toml["MATCHING"]["mz_tolerance"],
        "rt_tolerance": parsed_toml["MATCHING"]["rt_tolerance"],
        "rt_window": parsed_toml["MATCHING"]["ms1_rt_window"],
        "threshold": parsed_toml["FILTERING"]["total_intensity_threshold"],
    }


def __get_tmt_intensities(spectrum: Dict[str, Any]) -> Dict[str, float]:
    tmt_quants = {key: 0.0 for key in TMT.keys()}
    return tmt_quants


def __get_key(mass: float) -> int:
    return int(round(mass * 10000))


def __read_spectra(filename: str) -> Dict[str, Any]:
    spectra_ms1 = dict()
    spectra_ms2 = dict()
    total = 0
    total_ms1 = 0
    total_ms2 = 0
    with mzml.read(filename) as reader:
        for spectrum in reader:
            ms_level = int(spectrum["ms level"])
            if (
                "scanList" not in spectrum
                or "scan" not in spectrum["scanList"]
                or len(spectrum["scanList"]["scan"]) != 1
            ):
                raise RuntimeError(f"Can't get retention time for spectrum: {spectrum}")
            rt_in_min = float(spectrum["scanList"]["scan"][0]["scan start time"])
            rt_in_sec = rt_in_min * 60.0
            if ms_level == 2:
                if "precursorList" in spectrum:
                    for precursor in spectrum["precursorList"]["precursor"]:
                        for ion in precursor["selectedIonList"]["selectedIon"]:
                            primary_key = __get_key(float(ion["selected ion m/z"]))
                            secondary_key = __get_key(rt_in_sec)
                            s = dict()
                            s["precursor"] = float(ion["selected ion m/z"])
                            s["rt"] = rt_in_sec
                            s["mz_array"] = spectrum["m/z array"]
                            s["intensity_array"] = spectrum["intensity array"]
                            if primary_key in spectra_ms2:
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
                    raise KeyError(f"No precursor for MS2 spectrum found: {spectrum}")
            elif ms_level == 1:
                primary_key = __get_key(rt_in_sec)
                s = dict()
                s["precursor"] = None
                s["rt"] = rt_in_sec
                s["mz_array"] = spectrum["m/z array"]
                s["intensity_array"] = spectrum["intensity array"]
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


def __check_mz_in_ms1(
    precursor_mz: float, spectrum: Dict[str, Any], mz_tol: float
) -> bool:
    # todo
    return False


def __get_ms2_spectrum(
    precursor_mz: float | int,
    mz_tol: float,
    retention_time: float,
    rt_tol: float,
    retention_time_ms1_window: float,
    window_size_unidirectional: float,
    filter_threshold: float,
    spectra: Dict[str, Any],
) -> Dict[str, Any] | None:
    # todo
    primary_key_base = __get_key(precursor_mz)
    primary_key_window = __get_key(window_size_unidirectional)
    precursor = None
    for i in range(primary_key_window):
        if primary_key_base - i in spectra["ms2"]:
            precursor = spectra["ms2"][primary_key_base - i]
            break
        if primary_key_base + i in spectra["ms2"]:
            precursor = spectra["ms2"][primary_key_base + i]
            break
    if precursor is None:
        raise RuntimeError(
            f"Could not find a suitable precursor for precursor m/z {precursor_mz} and retention time {retention_time}."
        )
    secondary_key_base = __get_key(retention_time)
    secondary_key_window = __get_key(rt_tol)
    spectrum = None
    for i in range(secondary_key_window):
        if secondary_key_base - i in precursor:
            spectrum = precursor[secondary_key_base - i]
            break
        if secondary_key_base + i in precursor:
            spectrum = precursor[secondary_key_base + i]
            break
    if spectrum is None:
        raise RuntimeError(
            f"Could not find a suitable retention time for precursor m/z {precursor_mz} and retention time {retention_time}."
        )
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
    if ms1 is None:
        raise RuntimeError(
            f"Could not find a suitable MS1 spectrum for precursor m/z {precursor_mz} and retention time {retention_time}."
        )
    # intensity filter
    return


def __annotate_spectronaut_result(
    spectronaut_filename: str, spectra: Dict[str, Any], settings: Dict[str, Any]
) -> pd.DataFrame:
    df = pd.read_csv(spectronaut_filename, low_memory=False)
    channels = {key: [] for key in TMT.keys()}
    nr_of_filtered_spectra = 0
    for i, row in tqdm(
        df.iterrows(), total=df.shape[0], desc="Annotating Spectronaut result..."
    ):
        prec_mz = float(row[settings["precursor_mz"]])
        mz_tol = float(settings["mz_tolerance"])
        rt = float(row[settings["retention_time"]])
        rt_tol = float(settings["rt_tolerance"])
        window_size_unidirectional = float(settings["window_size"]) / 2.0
        filter_threshold = float(settings["threshold"])
        if not settings["retention_time_in_sec"]:
            rt = rt * 60.0
        spectrum = __get_ms2_spectrum(
            prec_mz,
            mz_tol,
            rt,
            rt_tol,
            window_size_unidirectional,
            filter_threshold,
            spectra,
        )
        if spectrum is not None:
            tmt_quants = __get_tmt_intensities(spectrum)
            for key in channels.keys():
                channels[key].append(tmt_quants[key])
        else:
            for key in channels.keys():
                channels[key].append(None)
            nr_of_filtered_spectra += 1
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
    parser.add_argument("--version", action="version", version=__version)
    args = parser.parse_args(argv)
    settings = __read_settings(args.config)
    spectra = __read_spectra(args.spectra)
    df = __annotate_spectronaut_result(args.spectronaut, spectra, settings)
    df.to_csv(args.spectronaut + "_tmt_quant.csv", index=False)
    return df


if __name__ == "__main__":
    _ = main()
