#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "pandas",
#   "numpy",
#   "tqdm",
#   "pyteomics[XML]",
#   "pyopenms",
# ]
# ///

# DIA TMT QUANTIFICATION CHIMERYS [MULTI-FILE]
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import os
import glob
from tmt_chimerys import __get_resolution_gui_map
from tmt_chimerys import __read_settings
from tmt_chimerys import __convert
from tmt_chimerys import __read_spectra_by_scannumber
from tmt_chimerys import __get_consensusXML_df
from tmt_chimerys import __get_consensusXML_map
from tmt_chimerys import __annotate_chimerys_result
from tmt_chimerys import __annotate_chimerys_protein_table
from tmt_chimerys import __annotate_result_conditions

CONFIG_FILE = "config.toml"
RESOLUTION_FILE = "resolution.csv"
WINDOW_FILE = None


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
        # run tmt chimerys
        settings = __read_settings(CONFIG_FILE)
        settings["window_size"] = float(w)
        print("Used settings:")
        print(settings)
        if WINDOW_FILE is not None:
            print(f"Using windows from given windows file: {WINDOW_FILE}")
        args_spectra = __convert(f"{f}.mzML")
        spectra = __read_spectra_by_scannumber(args_spectra)
        quantification_method = int(settings["quantification_method"])
        consensusXML_map = None
        if quantification_method != 1 and quantification_method != 3:
            consensusXML_df = __get_consensusXML_df(args_spectra)
            consensusXML_map = __get_consensusXML_map(consensusXML_df)
        df = __annotate_chimerys_result(
            filename=f"{f}_PSMs.txt",
            spectrum_filename=args_spectra,
            spectra=spectra,
            settings=settings,
            consensusXML_map=consensusXML_map,
            resolution_gui_map=resolution_gui_map,
            window_file=WINDOW_FILE,
        )
        df.to_csv(
            f"{f}_PSMs_purity_tmt_quant.txt",
            sep="\t",
            index=False,
        )
        df = __annotate_result_conditions(df, settings["conditions"])
        df.to_csv(
            f"{f}_PSMs_purity_tmt_quant_conditions.txt",
            sep="\t",
            index=False,
        )
        proteins_df = __annotate_chimerys_protein_table(
            f"{f}_Proteins.txt", df, settings
        )
        proteins_df.to_csv(
            f"{f}_Proteins_purity_tmt_quant.txt",
            sep="\t",
            index=False,
        )
        # console log
        print("---------- FINISHED ANALYSIS FOR ONE FILE ----------")
        print(f"processed {nr + 1}/{len(file_prefixes)} files.")


if __name__ == "__main__":
    main()
