#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "pandas",
#   "tqdm",
#   "pyteomics[XML]",
#   "pyopenms",
# ]
# ///

# DDA TMT QUANTIFICATION CHIMERYS [MULTI-FILE]
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import os
import glob
from tmt_chimerys import __get_resolution_gui_map
from tmt_chimerys import __read_settings
from tmt_chimerys import __convert
from tmt_spectronaut import __read_spectra
from tmt_chimerys import __get_consensusXML_df
from tmt_chimerys import __get_consensusXML_map
from tmt_diann import __annotate_diann_result

CONFIG_FILE = "config.toml"
RESOLUTION_FILE = "resolution.csv"
USE_OPENMS = True
# Verbose level where 0: ignore all warnings and errors, 1: raise warnings, and >= 2: raise errors
VERBOSE = 2


def main():
    mzml_files = [f for f in glob.glob("*.mzML")]
    file_prefixes = [os.path.splitext(f)[0] for f in mzml_files]
    print("Reading resolution.csv for all files...")
    resolution_gui_map = __get_resolution_gui_map(RESOLUTION_FILE)
    print("Sucessfully read resolution.csv for all files!")
    for nr, f in enumerate(file_prefixes):
        # try parsing window size just to check that it's correctly parsed
        # adjust if needed
        _ = float(f.split("DIA_mz")[1].split("_")[0].replace("c", "."))
        w = f.split("DIA_mz")[1].split("_")[0].replace("c", ".")
        # console log
        print("---------- STARTING ANALYSIS FOR ONE FILE ----------")
        print(f"file: {f}")
        print(f"window: {w}")
        # run tmt chimerys dda
        settings = __read_settings(CONFIG_FILE)
        settings["window_size"] = float(w)
        print("Used settings:")
        print(settings)
        args_spectra = __convert(f"{f}.mzML")
        spectra = __read_spectra(args_spectra)
        consensusXML_map = None
        if USE_OPENMS:
            consensusXML_df = __get_consensusXML_df(args_spectra)
            consensusXML_map = __get_consensusXML_map(consensusXML_df)
        df = __annotate_diann_result(
            diann_filename=f"{f}_precursor.parquet",
            spectrum_filename=args_spectra,
            spectra=spectra,
            settings=settings,
            consensusXML_map=consensusXML_map,
            resolution_gui_map=resolution_gui_map,
            verbose=VERBOSE,
        )
        df.to_parquet(
            f"{f}_precursor_purity_tmt_quant.parquet",
            index=False,
        )
        # console log
        print("---------- FINISHED ANALYSIS FOR ONE FILE ----------")
        print(f"processed {nr + 1}/{len(file_prefixes)} files.")


if __name__ == "__main__":
    main()
