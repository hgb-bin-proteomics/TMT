#!/usr/bin/env python3

# DIA TMT QUANTIFICATION DIA-NN
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import argparse
import pandas as pd
from tqdm import tqdm

from typing import Dict
from typing import Any

from tmt_chimerys import TMT
from tmt_chimerys import __get_bool_from_value
from tmt_chimerys import __read_settings
from tmt_chimerys import __get_tmt_intensities
from tmt_chimerys import __get_windows
from tmt_spectronaut import __read_spectra
from tmt_spectronaut import __get_ms2_spectrum

__version = "0.0.10"
__date = "2025-07-22"


# annotates the DIA-NN result with TMT quantities
# currently based on the report parquet
def __annotate_diann_result(
    diann_filename: str,
    spectra: Dict[str, Any],
    settings: Dict[str, Any],
    verbose: int = 2,
) -> pd.DataFrame:
    # spectra should be given by __read_spectra
    # settings should be given by __read_settings
    df = pd.read_parquet(diann_filename)
    channels = {key: [] for key in TMT.keys()}
    purities = list()
    nr_of_missing_ms1 = 0
    nr_of_impure_ids = 0
    for i, row in tqdm(
        df.iterrows(), total=df.shape[0], desc="Annotating DIA-NN result..."
    ):
        # get DIA-NN precursor m/z from identification
        prec_mz = float(row["Precursor.Mz"])
        # settings defined m/z tolerance
        mz_tol = float(settings["mz_tolerance"])
        # get DIA-NN retention time from identification
        rt = float(row["RT"]) * 60.0
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
        # isotope parameters
        do_deisotope = __get_bool_from_value(settings["deisotope"])
        isotope_tolerance = float(settings["isotope_tolerance"])
        max_charge = int(settings["max_charge"])
        # get corresponding MS2 spectrum for identification
        spectrum_purity = __get_ms2_spectrum(
            prec_mz,
            mz_tol,
            rt,
            rt_tol,
            rt_window,
            do_deisotope,
            isotope_tolerance,
            max_charge,
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
    # update DIA-NN result
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
        prog="tmt_diann.py",
        description="Calculates co-isolation purity for DIA-NN DIA TMT peptide matches and quantifies them.",
        epilog="(c) Research Institute of Molecular Pathology, 2025",
    )
    parser.add_argument(
        "-i",
        "--diann",
        dest="diann",
        required=True,
        help="Path/name of the DIA-NN result file.",
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
    df = __annotate_diann_result(args.diann, spectra, settings, int(args.verbose))
    df.to_parquet(args.diann.split(".parquet")[0] + "_purity_tmt_quant.parquet")
    return df


if __name__ == "__main__":
    _ = main()
