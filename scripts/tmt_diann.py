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
# ]
# ///

# DIA TMT QUANTIFICATION DIA-NN
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import argparse
import pandas as pd
from tqdm import tqdm

from typing import Optional
from typing import Dict
from typing import Any

from tmt_chimerys import TMT
from tmt_chimerys import RESOLUTION_GUI_COLS
from tmt_chimerys import __get_bool_from_value
from tmt_chimerys import __read_settings
from tmt_chimerys import __get_tmt_intensities
from tmt_chimerys import __get_tmt_intensities_oms
from tmt_chimerys import __get_consensusXML_df
from tmt_chimerys import __get_consensusXML_map
from tmt_chimerys import __get_resolution_gui_map
from tmt_chimerys import __get_resolution_gui_values
from tmt_chimerys import __get_tmt_intensities_resgui
from tmt_chimerys import __annotate_result_conditions
from tmt_chimerys import __subtract_noise
from tmt_chimerys import __get_windows
from tmt_chimerys import __convert
from tmt_spectronaut import __read_spectra
from tmt_spectronaut import __get_ms2_spectrum

__version = "2.0.1"
__date = "2025-11-19"


def __remove_ambiguous_pg(protein_table: pd.DataFrame) -> pd.DataFrame:
    protein_table["Filter:Is_Ambiguous_PG"] = protein_table.apply(
        lambda row: ";" in str(row["Protein.Group"]), axis=1
    )
    filtered_protein_table = protein_table[~protein_table["Filter:Is_Ambiguous_PG"]]
    if not isinstance(filtered_protein_table, pd.DataFrame):
        raise RuntimeError("Filtering did not return a table, too strict?")
    return filtered_protein_table


def __annotate_diann_pgs(
    precursor_table: pd.DataFrame,
    settings: Dict[str, Any],
) -> pd.DataFrame:
    q_value = float(settings["q_value"])
    min_reporter_res = float(settings["min_reporter_res"])
    min_purity = float(settings["min_purity"])
    conditions = settings["conditions"]
    has_resolution = "RESGUI_Resolution" in precursor_table.columns.tolist()
    psms_by_proteins = dict()
    for i, psm in tqdm(
        precursor_table.iterrows(),
        total=precursor_table.shape[0],
        desc="Filtering precursors...",
    ):
        pg = str(psm["Protein.Group"]).strip()
        purity = float(psm["Co-Isolation Purity"])
        global_q_value = float(psm["Global.Q.Value"])
        # remove PSMs below purity threshold
        if pd.isna(purity) or purity < min_purity:
            continue
        # remove PSMs above q-value
        if pd.isna(global_q_value) or global_q_value > q_value:
            continue
        if pg in psms_by_proteins:
            psms_by_proteins[pg].append(psm)
        else:
            psms_by_proteins[pg] = [psm]
    channels = {key: [] for key in TMT.keys()}
    for i, protein in tqdm(
        precursor_table.iterrows(),
        total=precursor_table.shape[0],
        desc="Annotating protein abundances...",
    ):
        pg = str(protein["Protein.Group"]).strip()
        psms_for_pg = list()
        if pg in psms_by_proteins:
            psms_for_pg = psms_by_proteins[pg]
        tmt_quants = {key: 0.0 for key in TMT.keys()}
        for psm in psms_for_pg:
            if has_resolution:
                for c in TMT.keys():
                    resgui_key = "RESGUI_" + c.split("-")[1] + " Resolution"
                    eligible_for_quant = True
                    if pd.isna(float(psm[resgui_key])):
                        eligible_for_quant = False
                    if float(psm[resgui_key]) < min_reporter_res:
                        eligible_for_quant = False
                    for condition in conditions:
                        if (
                            c in condition["reporters"]
                            and psm[f"Condition_SN_{condition['name']}"]
                            < condition["sn"]
                        ):
                            eligible_for_quant = False
                            break
                    if not eligible_for_quant:
                        tmt_quants[c] += 0.0
                    else:
                        tmt_quants[c] += psm[f"Annotated {c}"]
            else:
                for c in TMT.keys():
                    tmt_quants[c] += psm[f"Annotated {c}"]
        for k, v in tmt_quants.items():
            channels[k].append(v)
    for key in channels.keys():
        precursor_table[f"Annotated protein-level {key}"] = channels[key]
    return precursor_table


# annotates the DIA-NN result with TMT quantities
# currently based on the report parquet
def __annotate_diann_result(
    diann_filename: str,
    spectrum_filename: str,
    spectra: Dict[str, Any],
    settings: Dict[str, Any],
    consensusXML_map: Optional[Dict[int, Dict[int, pd.Series]]] = None,
    resolution_gui_map: Optional[Dict[str, Dict[int, pd.Series]]] = None,
    window_file: Optional[str] = None,
    verbose: int = 2,
) -> pd.DataFrame:
    # spectra should be given by __read_spectra
    # settings should be given by __read_settings
    df = pd.read_parquet(diann_filename)
    # subset to only precursors from ms file
    df = df[df["Run"] == spectrum_filename[:-5]]
    if not isinstance(df, pd.DataFrame) or df.shape[0] == 0:
        raise RuntimeError("Filtering for given MS file did not return a dataframe!")
    channels = {key: [] for key in TMT.keys()}
    resolution = {f"RESGUI_{key}": [] for key in RESOLUTION_GUI_COLS}
    purities = list()
    parsed_scannr = list()
    nr_of_missing_ms1 = 0
    nr_of_impure_ids = 0
    quantification_method = int(settings["quantification_method"])
    if quantification_method == 1:
        print("Using native quantification!")
    elif quantification_method == 3:
        print("Using Resolution GUI quantification!")
    else:
        print("Using OpenMS quantification!")
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
            float(settings["window_overlap"]),
            window_file,
        )
        # isotope parameters
        do_deisotope = __get_bool_from_value(settings["deisotope"])
        isotope_tolerance = float(settings["isotope_tolerance"])
        max_charge = int(settings["max_charge"])
        # normalization
        subtract_noise = __get_bool_from_value(settings["subtract_noise"])
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
            if quantification_method == 1:
                tmt_quants = __get_tmt_intensities(spectrum)
            elif quantification_method == 3:
                if resolution_gui_map is None:
                    raise ValueError(
                        "Quantification method Resolution GUI was selected but no resolution file was found!"
                    )
                else:
                    tmt_quants = __get_tmt_intensities_resgui(
                        spectrum, resolution_gui_map, spectrum_filename
                    )
            else:
                if consensusXML_map is None:
                    raise ValueError(
                        "Quantification method OpenMS was selected but no consensusXML was found!"
                    )
                else:
                    tmt_quants = __get_tmt_intensities_oms(spectrum, consensusXML_map)
            if subtract_noise:
                if resolution_gui_map is None:
                    raise ValueError(
                        "Subtract noise was set to true but no resolution file was found!"
                    )
                else:
                    tmt_quants = __subtract_noise(
                        tmt_quants, spectrum, resolution_gui_map, spectrum_filename
                    )
            for key in channels.keys():
                channels[key].append(tmt_quants[key])
            if resolution_gui_map is None:
                for key in resolution.keys():
                    resolution[key].append(None)
            else:
                resolution_values = __get_resolution_gui_values(
                    spectrum, resolution_gui_map, spectrum_filename
                )
                for key in resolution.keys():
                    resolution[key].append(resolution_values[key])
            purities.append(purity)
            if purity is None:
                nr_of_missing_ms1 += 1
            if purity is not None and purity < filter_threshold:
                nr_of_impure_ids += 1
            parsed_scannr.append(spectrum["scan_nr"])
        else:
            for key in channels.keys():
                channels[key].append(None)
            for key in resolution.keys():
                resolution[key].append(None)
            purities.append(None)
            nr_of_missing_ms1 += 1
            parsed_scannr.append(None)
    # update DIA-NN result
    df["Co-Isolation Purity"] = purities
    df["Parsed MS2 Scan Number"] = parsed_scannr
    for key in channels.keys():
        df[f"Annotated {key}"] = channels[key]
    if resolution_gui_map is not None:
        for key in resolution:
            df[key] = resolution[key]
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
        required=True,
        help="Path/name of the config file.",
        type=str,
    )
    parser.add_argument(
        "-r",
        "--resolution",
        dest="resolution",
        required=False,
        default=None,
        help="Path/name of the resolution.csv file from the Resolution GUI file.",
        type=str,
    )
    parser.add_argument(
        "-w",
        "--window",
        dest="window_file",
        default=None,
        help="Window file, overrides config file!",
        type=str,
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
    print(settings)
    if args.window_file is not None:
        print(f"Using windows from given windows file: {args.window_file}")
    args_spectra = __convert(args.spectra)
    spectra = __read_spectra(args_spectra)
    quantification_method = int(settings["quantification_method"])
    consensusXML_map = None
    if quantification_method != 1 and quantification_method != 3:
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
        window_file=args.window_file,
        verbose=int(args.verbose),
    )
    df = __annotate_result_conditions(df, settings["conditions"])
    df = __annotate_diann_pgs(df, settings)
    if not __get_bool_from_value(settings["keep_pg"]):
        df = __remove_ambiguous_pg(df)
    df.to_parquet(
        args.proteins.split(".parquet")[0] + "_purity_tmt_quant.parquet",
        index=False,
    )
    return df


if __name__ == "__main__":
    _ = main()
