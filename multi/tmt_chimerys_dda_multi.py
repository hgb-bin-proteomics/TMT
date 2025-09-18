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
from tmt_chimerys_dda import main as tmt_chimerys_dda

RESOLUTION_FILE = "resolution.csv"


def main():
    mzml_files = [f for f in glob.glob("*.mzML")]
    file_prefixes = [os.path.splitext(f)[0] for f in mzml_files]
    for nr, f in enumerate(file_prefixes):
        # console log
        print("---------- STARTING ANALYSIS FOR ONE FILE ----------")
        print(f"file: {f}")
        # call tmt chimerys script
        _ = tmt_chimerys_dda(
            [
                "-i",
                f"{f}_PSMs.txt",
                "-s",
                f"{f}.mzML",
                "-c",
                "config.toml",
                "-p",
                f"{f}_Proteins.txt",
                "-r",
                RESOLUTION_FILE,
            ]
        )
        # console log
        print("---------- FINISHED ANALYSIS FOR ONE FILE ----------")
        print(f"processed {nr + 1}/{len(file_prefixes)} files.")


if __name__ == "__main__":
    main()
