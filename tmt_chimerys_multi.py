#!/usr/bin/env python3

# UNNAMED TMT PROJECT
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import os
import glob
from tmt_chimerys import main as tmt_chimerys


def main():
    mzml_files = [f for f in glob.glob("*.mzML")]
    file_prefixes = [os.path.splitext(f)[0] for f in mzml_files]
    for f in file_prefixes:
        # try parsing window size just to check that it's correctly parsed
        _ = float(f.split("pool_DIA_mz")[1].split("_")[0].replace("c", "."))
        w = f.split("pool_DIA_mz")[1].split("_")[0].replace("c", ".")
        # call tmt chimerys script
        _ = tmt_chimerys(
            ["-i", f"{f}_PSMs.txt", "-s", f"{f}.mzML", "-c", "config.toml", "-w", w]
        )
