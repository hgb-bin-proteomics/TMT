#!/usr/bin/env python3

# UNNAMED TMT PROJECT
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

__version = "0.0.1"
__date = "2025-06-25"

import tomllib
import argparse
import pandas as pd


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
    return mz * charge - proton * charge


def __get_settings(toml: str) -> dict:
    parsed_toml = None
    with open(toml, "rb") as f:
        parsed_toml = tomllib.load(f)
        f.close()
    if parsed_toml is None:
        raise RuntimeError()
    return {
        "precursor_mass": parsed_toml["SPECTRONAUT"]["precursor_mass"],
        "precursor_mz": parsed_toml["SPECTRONAUT"]["precursor_mz"],
        "precursor_charge": parsed_toml["SPECTRONAUT"]["precursor_mz"],
        "tolerance": parsed_toml["MATCHING"]["tolerance"],
    }


def __get_tmt_intensities(spectrum: dict) -> dict:
    return


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(
        prog="ProgramName",
        description="What the program does",
        epilog="Text at the bottom of help",
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
    settings = __get_settings(args.config)
    return
