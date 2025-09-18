#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "pandas",
#   "pyarrow",
#   "tqdm",
#   "pyteomics[XML]",
#   "pyopenms",
#   "gooey",
# ]
# ///

# DIA TMT QUANTIFICATION DIA-NN
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pandas as pd

from gooey import Gooey
from gooey import GooeyParser

from tmt_chimerys import __read_settings
from tmt_chimerys import __get_consensusXML_df
from tmt_chimerys import __get_consensusXML_map
from tmt_chimerys import __get_resolution_gui_map
from tmt_chimerys import __convert
from tmt_spectronaut import __read_spectra
from tmt_diann import __annotate_diann_result

__version = "1.0.2"
__date = "2025-09-18"


@Gooey(
    encoding="utf-8",
    program_name=f"TMT DIA-NN {__version}",
    menu=[
        {
            "name": "Help",
            "items": [
                {
                    "type": "Link",
                    "menuTitle": "Project Page",
                    "url": "https://github.com/hgb-bin-proteomics/TMT/",
                }
            ],
        }
    ],
)
def main(argv=None) -> pd.DataFrame:
    parser = GooeyParser(
        prog="tmt_diann.py",
        description="Calculates co-isolation purity for DIA-NN DIA TMT peptide matches and quantifies them.",
        epilog="(c) Research Institute of Molecular Pathology, 2025",
    )
    req = parser.add_argument_group("Required", "Required Arguments.")
    req.add_argument(
        "-i",
        "--diann",
        dest="diann",
        required=True,
        help="Path/name of the DIA-NN result file.",
        type=str,
        widget="FileChooser",
    )
    req.add_argument(
        "-s",
        "--spectra",
        dest="spectra",
        required=True,
        help="Path/name of the mass spectra file in mzML format.",
        type=str,
        widget="FileChooser",
    )
    req.add_argument(
        "-c",
        "--config",
        dest="config",
        default=None,
        help="Path/name of the config file.",
        type=str,
        widget="FileChooser",
    )
    opt = parser.add_argument_group("Optional", "Optional Arguments.")
    opt.add_argument(
        "-r",
        "--resolution",
        dest="resolution",
        required=False,
        default=None,
        help="Path/name of the resolution.csv file from the Resolution GUI file.",
        type=str,
        widget="FileChooser",
    )
    opt.add_argument(
        "-w",
        "--window",
        dest="window",
        default=None,
        help="Window size, overrides config file!",
        type=float,
    )
    opt.add_argument(
        "-n",
        "--native",
        dest="native",
        default=False,
        action="store_true",
        help="Use native quantification instead of OpenMS quantification which is used by default.",
    )
    opt.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        default=2,
        help="Verbose level.",
        type=int,
        widget="IntegerField",
        gooey_options={"initial_value": 1, "min": 0, "max": 2, "increment": 1},
    )
    args = parser.parse_args(argv)
    settings = __read_settings(args.config)
    if args.window is not None:
        settings["window_size"] = float(args.window)
    print("Read settings:")
    print(settings)
    args_spectra = __convert(args.spectra)
    spectra = __read_spectra(args_spectra)
    consensusXML_map = None
    if not args.native:
        consensusXML_df = __get_consensusXML_df(args_spectra)
        consensusXML_map = __get_consensusXML_map(consensusXML_df)
    resolution_gui_map = None
    if args.resolution is not None:
        resolution_gui_map = __get_resolution_gui_map(args.resolution)
    df = __annotate_diann_result(
        diann_filename=args.diann,
        spectrum_filename=args_spectra,
        spectra=spectra,
        settings=settings,
        consensusXML_map=consensusXML_map,
        resolution_gui_map=resolution_gui_map,
        verbose=int(args.verbose),
    )
    df.to_parquet(
        args.diann.split(".parquet")[0] + "_purity_tmt_quant.parquet", index=False
    )
    print("Script finished successfully!")
    return df


if __name__ == "__main__":
    _ = main()
