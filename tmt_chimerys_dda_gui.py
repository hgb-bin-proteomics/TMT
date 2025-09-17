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

# DDA TMT QUANTIFICATION CHIMERYS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pandas as pd

from gooey import Gooey
from gooey import GooeyParser

from tmt_chimerys import __read_settings
from tmt_chimerys import __read_spectra_by_scannumber
from tmt_chimerys import __get_consensusXML_df
from tmt_chimerys import __get_consensusXML_map
from tmt_chimerys import __get_resolution_gui_map
from tmt_chimerys import __annotate_chimerys_protein_table
from tmt_chimerys import __convert
from tmt_chimerys_dda import __annotate_chimerys_result

__version = "1.0.0"
__date = "2025-09-17"


@Gooey(
    encoding="utf-8",
    program_name="TMT Chimerys",
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
        prog="tmt_chimerys_dda.py",
        description="Calculates co-isolation purity for Chimerys DDA TMT PSMs and optionally quantifies them.",
        epilog="(c) Research Institute of Molecular Pathology, 2025",
    )
    req = parser.add_argument_group("Required", "Required Arguments.")
    req.add_argument(
        "-i",
        "--chimerys",
        dest="chimerys",
        required=True,
        help="Path/name of the Chimerys result file in tab-separated .txt format.",
        type=str,
    )
    req.add_argument(
        "-s",
        "--spectra",
        dest="spectra",
        required=True,
        help="Path/name of the mass spectra file in mzML format.",
        type=str,
    )
    req.add_argument(
        "-c",
        "--config",
        dest="config",
        required=True,
        help="Path/name of the config file.",
        type=str,
    )
    opt = parser.add_argument_group("Optional", "Optional Arguments.")
    opt.add_argument(
        "-p",
        "--proteins",
        dest="proteins",
        required=False,
        default=None,
        help="Path/name of the Chimerys protein result file in tab-separated .txt format.",
        type=str,
    )
    opt.add_argument(
        "-r",
        "--resolution",
        dest="resolution",
        required=False,
        default=None,
        help="Path/name of the resolution.csv file from the Resolution GUI file.",
        type=str,
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
    args = parser.parse_args(argv)
    settings = __read_settings(args.config)
    if args.window is not None:
        settings["window_size"] = float(args.window)
    print("Read settings:")
    print(settings)
    args_spectra = __convert(args.spectra)
    spectra = __read_spectra_by_scannumber(args_spectra)
    consensusXML_map = None
    if not args.native:
        consensusXML_df = __get_consensusXML_df(args_spectra)
        consensusXML_map = __get_consensusXML_map(consensusXML_df)
    resolution_gui_map = None
    if args.resolution is not None:
        resolution_gui_map = __get_resolution_gui_map(args.resolution)
    df = __annotate_chimerys_result(
        filename=args.chimerys,
        spectrum_filename=args_spectra,
        spectra=spectra,
        settings=settings,
        consensusXML_map=consensusXML_map,
        resolution_gui_map=resolution_gui_map,
    )
    df.to_csv(
        args.chimerys.split(".txt")[0] + "_purity_tmt_quant.txt",
        sep="\t",
        index=False,
    )
    if args.proteins is not None:
        proteins_df = __annotate_chimerys_protein_table(args.proteins, df, settings)
        proteins_df.to_csv(
            args.proteins.split(".txt")[0] + "_purity_tmt_quant.txt",
            sep="\t",
            index=False,
        )
    print("Script finished successfully!")
    return df


if __name__ == "__main__":
    _ = main()
