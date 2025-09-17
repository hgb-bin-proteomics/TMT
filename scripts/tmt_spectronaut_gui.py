#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "pandas",
#   "tqdm",
#   "pyteomics[XML]",
#   "pyopenms",
#   "gooey",
# ]
# ///

# DIA TMT QUANTIFICATION SPECTRONAUT
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
from tmt_spectronaut import __annotate_spectronaut_result

__version = "1.0.1"
__date = "2025-09-17"


@Gooey(
    encoding="utf-8",
    program_name="TMT Spectronaut",
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
        prog="tmt_spectronaut.py",
        description="Calculates co-isolation purity for Spectronaut DIA TMT peptide matches and quantifies them.",
        epilog="(c) Research Institute of Molecular Pathology, 2025",
    )
    req = parser.add_argument_group("Required", "Required Arguments.")
    req.add_argument(
        "-i",
        "--spectronaut",
        dest="spectronaut",
        required=True,
        help="Path/name of the Spectronaut result file.",
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
        help="Verbose level where 0: ignore all warnings and errors, 1: raise warnings, and >= 2: raise errors!",
        type=int,
        widget="IntegerField",
        gooey_options={"initial_value": 1, "min": 0, "max": 2, "increment": 1},
    )
    args = parser.parse_args(argv)
    settings = __read_settings(args.config)
    args_spectra = __convert(args.spectra)
    if args.window is not None:
        settings["window_size"] = float(args.window)
    print("Read settings:")
    print(settings)
    spectra = __read_spectra(args_spectra)
    consensusXML_map = None
    if not args.native:
        consensusXML_df = __get_consensusXML_df(args_spectra)
        consensusXML_map = __get_consensusXML_map(consensusXML_df)
    resolution_gui_map = None
    if args.resolution is not None:
        resolution_gui_map = __get_resolution_gui_map(args.resolution)
    df = __annotate_spectronaut_result(
        spectronaut_filename=args.spectronaut,
        spectrum_filename=args_spectra,
        spectra=spectra,
        settings=settings,
        consensusXML_map=consensusXML_map,
        resolution_gui_map=resolution_gui_map,
        verbose=int(args.verbose),
    )
    df.to_csv(
        args.spectronaut.split(".csv")[0] + "_purity_tmt_quant.csv",
        sep=",",
        index=False,
    )
    print("Script finished successfully!")
    return df


if __name__ == "__main__":
    _ = main()
