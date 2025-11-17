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

# DIA TMT QUANTIFICATION SPECTRONAUT
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import argparse
import warnings
import pandas as pd
from tqdm import tqdm
from pyteomics import mzml

from typing import Optional
from typing import Dict
from typing import List
from typing import Tuple
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
from tmt_chimerys import __get_key
from tmt_chimerys import __parse_scan_nr_from_id
from tmt_chimerys import __check_mz_in_ms1
from tmt_chimerys import __calculate_precursor_intensity_ms1
from tmt_chimerys import __get_windows
from tmt_chimerys import __convert

__version = "2.0.0"
__date = "2025-11-17"


def __remove_ambiguous_pg(protein_table: pd.DataFrame) -> pd.DataFrame:
    protein_table["Filter:Is_Ambiguous_PG"] = protein_table.apply(
        lambda row: ";" in str(row["PG.ProteinGroups"]), axis=1
    )
    filtered_protein_table = protein_table[~protein_table["Filter:Is_Ambiguous_PG"]]
    if not isinstance(filtered_protein_table, pd.DataFrame):
        raise RuntimeError("Filtering did not return a table, too strict?")
    return filtered_protein_table


def __annotate_spectronaut_protein_df(
    protein_table: pd.DataFrame,
    psm_table: pd.DataFrame,
    min_reporter_res: float,
    min_purity: float,
    conditions: List[Dict[str, Any]],
) -> pd.DataFrame:
    has_resolution = "RESGUI_Resolution" in psm_table.columns.tolist()
    psms_by_proteins = dict()
    for i, psm in tqdm(
        psm_table.iterrows(), total=psm_table.shape[0], desc="Filtering precursors..."
    ):
        pg = str(psm["PG.ProteinGroups"]).strip()
        purity = float(psm["Co-Isolation Purity"])
        # remove PSMs below purity threshold
        if pd.isna(purity) or purity < min_purity:
            continue
        if pg in psms_by_proteins:
            psms_by_proteins[pg].append(psm)
        else:
            psms_by_proteins[pg] = [psm]
    channels = {key: [] for key in TMT.keys()}
    for i, protein in tqdm(
        protein_table.iterrows(),
        total=protein_table.shape[0],
        desc="Annotating protein abundances...",
    ):
        pg = str(protein["PG.ProteinGroups"]).strip()
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
        protein_table[f"Annotated protein-level {key}"] = channels[key]
    return protein_table


def __annotate_spectronaut_protein_table(
    protein_table: str,
    psm_table: pd.DataFrame,
    settings: Dict[str, Any],
) -> pd.DataFrame:
    protein_df = pd.read_csv(protein_table, low_memory=False)
    min_reporter_res = float(settings["min_reporter_res"])
    min_purity = float(settings["min_purity"])
    return __annotate_spectronaut_protein_df(
        protein_df,
        psm_table,
        min_reporter_res=min_reporter_res,
        min_purity=min_purity,
        conditions=settings["conditions"],
    )


# read mass spectra from an mzML file
def __read_spectra(filename: str) -> Dict[str, Any]:
    # reads an mzML file and sorts spectra into an easily accessible data structure
    # that stores them based on precursor m/z and retention time
    # returned is a Dict[str, Dict[int, Any]]
    # at the top level there are two keys: "ms1" whichs maps to a dictionary containing
    # all MS1 spectra
    # and "ms2" which maps to a dictionary containing all MS2 spectra
    # MS1: Dict[int, Spectrum]
    # maps retention time in seconds * 10 000 (rounded) to MS1 spectra
    # MS2: Dict[int, Dict[int, Spectrum]]
    # maps precursor m/z * 10 000 (rounded) to dictionaries that map retention time
    # in seconds * 10 000 (rounded) to MS2 spectra
    spectra_ms1 = dict()
    spectra_ms2 = dict()
    total = 0
    total_ms1 = 0
    total_ms2 = 0
    with mzml.read(filename) as reader:
        for spectrum in reader:
            scan = __parse_scan_nr_from_id(spectrum["id"])
            # get MS level
            ms_level = int(spectrum["ms level"])
            # check if all fields for retrieving retention time are available
            if (
                "scanList" not in spectrum
                or "scan" not in spectrum["scanList"]
                or len(spectrum["scanList"]["scan"]) != 1
            ):
                raise RuntimeError(f"Can't get retention time for spectrum: {spectrum}")
            # get retention time
            rt_in_min = float(spectrum["scanList"]["scan"][0]["scan start time"])
            # retention time seems to be in minutes in mzML files
            rt_in_sec = rt_in_min * 60.0
            # for MS2 spectra we extract all precursors
            # though realistically there should only be one
            if ms_level == 2:
                if "precursorList" not in spectrum:
                    raise RuntimeError(
                        f"[precursorList] No precursor for MS2 spectrum found: {spectrum}"
                    )
                if (
                    "precursor" not in spectrum["precursorList"]
                    or len(spectrum["precursorList"]["precursor"]) != 1
                ):
                    raise RuntimeError(
                        f"[precursor] No precursor for MS2 spectrum found: {spectrum}"
                    )
                for precursor in spectrum["precursorList"]["precursor"]:
                    if "selectedIonList" not in precursor:
                        raise RuntimeError(
                            f"[selectedIonList] No precursor for MS2 spectrum found: {spectrum}"
                        )
                    if (
                        "selectedIon" not in precursor["selectedIonList"]
                        or len(precursor["selectedIonList"]["selectedIon"]) != 1
                    ):
                        raise RuntimeError(
                            f"[selectedIon] No precursor for MS2 spectrum found: {spectrum}"
                        )
                    for ion in precursor["selectedIonList"]["selectedIon"]:
                        # the primary key is the precursor m/z * 10 000 (rounded)
                        primary_key = __get_key(float(ion["selected ion m/z"]))
                        # the secondary key is the retention time in seconds * 10 000 (rounded)
                        secondary_key = __get_key(rt_in_sec)
                        s = dict()
                        s["scan_nr"] = scan
                        s["precursor"] = float(ion["selected ion m/z"])
                        s["rt"] = rt_in_sec
                        s["mz_array"] = spectrum["m/z array"]
                        s["intensity_array"] = spectrum["intensity array"]
                        if primary_key in spectra_ms2:
                            # if there already is another MS2 spectrum with the same precursor m/z
                            # and same retention time an error is raised
                            if secondary_key in spectra_ms2[primary_key]:
                                raise KeyError(
                                    f"MS2 spectrum for precursor {s['precursor']} and retention time {s['rt']} already exists!"
                                )
                            else:
                                spectra_ms2[primary_key][secondary_key] = s
                                total += 1
                                total_ms2 += 1
                        else:
                            spectra_ms2[primary_key] = {secondary_key: s}
                            total += 1
                            total_ms2 += 1
            # for MS1 spectra we save them based on retention time
            elif ms_level == 1:
                # the primary key is the retention time in seconds * 10 000 (rounded)
                primary_key = __get_key(rt_in_sec)
                s = dict()
                s["scan_nr"] = scan
                s["precursor"] = None
                s["rt"] = rt_in_sec
                s["mz_array"] = spectrum["m/z array"]
                s["intensity_array"] = spectrum["intensity array"]
                # if there are multiple MS1 spectra with the same retention time
                # an error is raised
                if primary_key in spectra_ms1:
                    raise KeyError(
                        f"MS1 spectrum for retention time {s['rt']} already exists!"
                    )
                else:
                    spectra_ms1[primary_key] = s
                    total += 1
                    total_ms1 += 1
            else:
                raise ValueError(
                    f"Found spectrum with MS level {ms_level}. Not supported!"
                )
    print(f"Total number of parsed spectra: {total}")
    print(f"Number of MS1 spectra: {total_ms1}")
    print(f"Number of MS2 spectra: {total_ms2}")
    return {"ms1": spectra_ms1, "ms2": spectra_ms2}


# get corresponding MS2 spectrum for a Spectronaut given precursor m/z and
# retention time
# May return None if no MS2 spectrum or corresponding MS1 spectrum is found
def __get_ms2_spectrum(
    precursor_mz: float,
    mz_tol: float,
    retention_time: float,
    rt_tol: float,
    retention_time_ms1_window: float,
    do_deisotope: bool,
    isotope_tol: float,
    max_charge: int,
    window_size_unidirectional: float,
    noise_threshold: float,
    spectra: Dict[str, Any],
    windows: List[Tuple[float, float]],
    verbose: int,
) -> Dict[str, Any]:
    # get by most similar precursor m/z within window
    primary_key_base = __get_key(precursor_mz)
    primary_key_window = __get_key(window_size_unidirectional)
    precursor = None
    # take the closest precursor m/z MS2 spectrum within the window_size (half)
    for i in range(primary_key_window):
        if primary_key_base - i in spectra["ms2"]:
            precursor = spectra["ms2"][primary_key_base - i]
            break
        if primary_key_base + i in spectra["ms2"]:
            precursor = spectra["ms2"][primary_key_base + i]
            break
    # if no precursor is found, raise an error
    if precursor is None:
        if verbose > 1:
            raise RuntimeError(
                f"[no precursor] Could not find a suitable precursor for precursor m/z {precursor_mz} and retention time {retention_time}."
            )
        elif verbose == 1:
            warnings.warn(
                RuntimeWarning(
                    f"[no precursor] Could not find a suitable precursor for precursor m/z {precursor_mz} and retention time {retention_time}."
                )
            )
            return {"spectrum": None, "purity": None}
        else:
            return {"spectrum": None, "purity": None}
    # get by most similar retention time using a retention time tolerance
    secondary_key_base = __get_key(retention_time)
    secondary_key_window = __get_key(rt_tol)
    spectrum = None
    # take the MS2 spectrum closest in retention time but maximum rt_tol seconds away
    for i in range(secondary_key_window):
        if secondary_key_base - i in precursor:
            spectrum = precursor[secondary_key_base - i]
            break
        if secondary_key_base + i in precursor:
            spectrum = precursor[secondary_key_base + i]
            break
    # if no MS2 spectrum is found an error is raised
    if spectrum is None:
        if verbose > 1:
            rt_deltas = list()
            for rt in precursor.keys():
                rt_deltas.append(abs(secondary_key_base - rt))
            rt_deltas = sorted(rt_deltas)
            raise RuntimeError(
                f"[no rt] Could not find a suitable retention time for precursor m/z {precursor_mz} and retention time {retention_time}."
                f" Closest retention time: {rt_deltas[0] / 10000.0} s"
            )
        elif verbose == 1:
            rt_deltas = list()
            for rt in precursor.keys():
                rt_deltas.append(abs(secondary_key_base - rt))
            rt_deltas = sorted(rt_deltas)
            warnings.warn(
                RuntimeWarning(
                    f"[no rt] Could not find a suitable retention time for precursor m/z {precursor_mz} and retention time {retention_time}."
                    f" Closest retention time: {rt_deltas[0] / 10000.0} s"
                )
            )
            return {"spectrum": None, "purity": None}
        else:
            return {"spectrum": None, "purity": None}
    # look for closest MS1 spectrum that has precursor (with m/z tolerance mz_tol)
    # within ms1 rt window
    retention_time_key = __get_key(spectrum["rt"])
    retention_time_ms1_window_range = __get_key(retention_time_ms1_window)
    ms1 = None
    for i in range(retention_time_ms1_window_range):
        if retention_time_key - i in spectra["ms1"]:
            if __check_mz_in_ms1(
                precursor_mz, spectra["ms1"][retention_time_key - i], mz_tol
            ):
                ms1 = spectra["ms1"][retention_time_key - i]
                break
        if retention_time_key + i in spectra["ms1"]:
            if __check_mz_in_ms1(
                precursor_mz, spectra["ms1"][retention_time_key + i], mz_tol
            ):
                ms1 = spectra["ms1"][retention_time_key + i]
                break
    # if no MS1 spectrum is found -> no error, but None returned
    if ms1 is None:
        warnings.warn(
            RuntimeWarning(
                f"Could not find a suitable MS1 spectrum for precursor m/z {precursor_mz} and retention time {spectrum['rt']}."
            )
        )
        return {"spectrum": spectrum, "purity": None}
    # intensity filter
    purity = __calculate_precursor_intensity_ms1(
        precursor_mz,
        ms1,
        spectrum,
        mz_tol,
        do_deisotope,
        isotope_tol,
        max_charge,
        noise_threshold,
        windows,
    )
    return {"spectrum": spectrum, "purity": purity}


# annotates the Spectronaut result with TMT quantities
# currently based on the factory report
def __annotate_spectronaut_result(
    spectronaut_filename: str,
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
    df = pd.read_csv(spectronaut_filename, low_memory=False)
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
        df.iterrows(), total=df.shape[0], desc="Annotating Spectronaut result..."
    ):
        # get Spectronaut precursor m/z from identification
        prec_mz = float(row["FG.PrecMz"])
        # settings defined m/z tolerance
        mz_tol = float(settings["mz_tolerance"])
        # get Spectronaut retention time from identification
        rt = float(row["EG.ApexRT"]) * 60.0
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
    # update Spectronaut result
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
        prog="tmt_spectronaut.py",
        description="Calculates co-isolation purity for Spectronaut DIA TMT peptide matches and quantifies them.",
        epilog="(c) Research Institute of Molecular Pathology, 2025",
    )
    parser.add_argument(
        "-i",
        "--spectronaut",
        dest="spectronaut",
        required=True,
        help="Path/name of the Spectronaut result file.",
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
        "-p",
        "--proteins",
        dest="proteins",
        required=False,
        default=None,
        help="Path/name of the Spectronaut protein result file in comma-separated .csv format.",
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
        help="Verbose level where 0: ignore all warnings and errors, 1: raise warnings, and >= 2: raise errors!",
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
    df = __annotate_spectronaut_result(
        spectronaut_filename=args.spectronaut,
        spectrum_filename=args_spectra,
        spectra=spectra,
        settings=settings,
        consensusXML_map=consensusXML_map,
        resolution_gui_map=resolution_gui_map,
        window_file=args.window_file,
        verbose=int(args.verbose),
    )
    df.to_csv(
        args.spectronaut.split(".csv")[0] + "_purity_tmt_quant.csv",
        sep=",",
        index=False,
    )
    df = __annotate_result_conditions(df, settings["conditions"])
    df.to_csv(
        args.spectronaut.split(".csv")[0] + "_purity_tmt_quant_conditions.csv",
        sep=",",
        index=False,
    )
    if args.proteins is not None:
        proteins_df = __annotate_spectronaut_protein_table(args.proteins, df, settings)
        if not __get_bool_from_value(settings["keep_pg"]):
            proteins_df = __remove_ambiguous_pg(proteins_df)
        proteins_df.to_csv(
            args.proteins.split(".csv")[0] + "_purity_tmt_quant_pg.csv",
            sep=",",
            index=False,
        )
    return df


if __name__ == "__main__":
    _ = main()
