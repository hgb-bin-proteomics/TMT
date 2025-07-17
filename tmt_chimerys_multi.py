#!/usr/bin/env python3

# DIA TMT QUANTIFICATION CHIMERYS [MULTI-FILE]
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import os
import glob
from tmt_chimerys import main as tmt_chimerys


def main():
    mzml_files = [f for f in glob.glob("*.mzML")]
    file_prefixes = [os.path.splitext(f)[0] for f in mzml_files]
    for nr, f in enumerate(file_prefixes):
        # try parsing window size just to check that it's correctly parsed
        _ = float(f.split("pool_DIA_mz")[1].split("_")[0].replace("c", "."))
        w = f.split("pool_DIA_mz")[1].split("_")[0].replace("c", ".")
        # console log
        print("---------- STARTING ANALYSIS FOR ONE FILE ----------")
        print(f"file: {f}")
        print(f"window: {w}")
        # call tmt chimerys script
        _ = tmt_chimerys(
            ["-i", f"{f}_PSMs.txt", "-s", f"{f}.mzML", "-c", "config.toml", "-w", w]
        )
        # console log
        print("---------- FINISHED ANALYSIS FOR ONE FILE ----------")
        print(f"processed {nr + 1}/{len(file_prefixes)} files.")


if __name__ == "__main__":
    main()
