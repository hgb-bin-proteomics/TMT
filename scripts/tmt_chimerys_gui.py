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

# DIA TMT QUANTIFICATION CHIMERYS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import os
import tomllib
import warnings
import subprocess
import pandas as pd
from tqdm import tqdm
from pyteomics import mzml
import pyopenms as oms
import urllib.request
import zipfile

from typing import Optional
from typing import Dict
from typing import List
from typing import Tuple
from typing import Any

from gooey import Gooey
from gooey import GooeyParser


__version = "2.1.1"
__date = "2025-11-20"

TMT_TOLERANCE = 0.0025
TMT = {
    "TMTpro-126": 126.127726,
    "TMTpro-127N": 127.124761,
    "TMTpro-127C": 127.131081,
    "TMTpro-128N": 128.128116,
    "TMTpro-128C": 128.134436,
    "TMTpro-129N": 129.131471,
    "TMTpro-129C": 129.13779,
    "TMTpro-130N": 130.134825,
    "TMTpro-130C": 130.141145,
    "TMTpro-131N": 131.13818,
    "TMTpro-131C": 131.1445,
    "TMTpro-132N": 132.141535,
    "TMTpro-132C": 132.147855,
    "TMTpro-133N": 133.14489,
    "TMTpro-133C": 133.15121,
    "TMTpro-134N": 134.148245,
    "TMTpro-134C": 134.154565,
    "TMTpro-135N": 135.151600,
}
TMT_OMS = {
    "TMTpro-126": "tmt18plex_126",
    "TMTpro-127N": "tmt18plex_127N",
    "TMTpro-127C": "tmt18plex_127C",
    "TMTpro-128N": "tmt18plex_128N",
    "TMTpro-128C": "tmt18plex_128C",
    "TMTpro-129N": "tmt18plex_129N",
    "TMTpro-129C": "tmt18plex_129C",
    "TMTpro-130N": "tmt18plex_130N",
    "TMTpro-130C": "tmt18plex_130C",
    "TMTpro-131N": "tmt18plex_131N",
    "TMTpro-131C": "tmt18plex_131C",
    "TMTpro-132N": "tmt18plex_132N",
    "TMTpro-132C": "tmt18plex_132C",
    "TMTpro-133N": "tmt18plex_133N",
    "TMTpro-133C": "tmt18plex_133C",
    "TMTpro-134N": "tmt18plex_134N",
    "TMTpro-134C": "tmt18plex_134C",
    "TMTpro-135N": "tmt18plex_135N",
}
RESOLUTION_GUI_COLS = [
    "Resolution",
    "TIC",
    "126 Resolution",
    "126 Intensity",
    "126 Noise",
    "127N Resolution",
    "127N Intensity",
    "127N Noise",
    "127C Resolution",
    "127C Intensity",
    "127C Noise",
    "128N Resolution",
    "128N Intensity",
    "128N Noise",
    "128C Resolution",
    "128C Intensity",
    "128C Noise",
    "129N Resolution",
    "129N Intensity",
    "129N Noise",
    "129C Resolution",
    "129C Intensity",
    "129C Noise",
    "130N Resolution",
    "130N Intensity",
    "130N Noise",
    "130C Resolution",
    "130C Intensity",
    "130C Noise",
    "131N Resolution",
    "131N Intensity",
    "131N Noise",
    "131C Resolution",
    "131C Intensity",
    "131C Noise",
    "132N Resolution",
    "132N Intensity",
    "132N Noise",
    "132C Resolution",
    "132C Intensity",
    "132C Noise",
    "133N Resolution",
    "133N Intensity",
    "133N Noise",
    "133C Resolution",
    "133C Intensity",
    "133C Noise",
    "134N Resolution",
    "134N Intensity",
    "134N Noise",
    "134C Resolution",
    "134C Intesntiy",
    "134C Noise",
    "135N Resolution",
    "135N Intensity",
    "135N Noise",
]
PROTON = 1.007276466812
ISOTOPE = 1.00335
STRATEGY = 1


def __convert(filename: str) -> str:
    dl_url = "https://github.com/CompOmics/ThermoRawFileParser/releases/download/v1.4.5/ThermoRawFileParser1.4.5.zip"
    if filename[-5:].lower() == ".mzml":
        print(f"Found mzML file with name {filename}.")
        print("Not converting file...")
        return filename
    print(f"Found RAW file with name {filename}.")
    print("Converting using ThermoRawFileParser...")
    if os.path.isdir("ThermoRawFileParser1.4.5") and os.path.exists(
        os.path.join("ThermoRawFileParser1.4.5", "ThermoRawFileParser.exe")
    ):
        print("Found existing ThermoRawFileParser installation!")
    else:
        print("Downloading ThermoRawFileParser!")
        urllib.request.urlretrieve(dl_url, "ThermoRawFileParser1.4.5.zip")
        with zipfile.ZipFile("ThermoRawFileParser1.4.5.zip", "r") as f:
            f.extractall("ThermoRawFileParser1.4.5")
    subprocess.call(
        [
            "ThermoRawFileParser1.4.5/ThermoRawFileParser.exe",
            "-i",
            filename,
        ]
    )
    new_filename = filename[:-3] + "mzML"
    print(f"Successfully converted file to {new_filename}!")
    return new_filename


def __get_sn_for_condition(row: pd.Series, reporters: List[str]) -> float:
    total_signal = 0.0
    total_noise = 0.0
    for reporter in reporters:
        reporter_signal = 0.0
        reporter_noise = 0.0
        label = reporter.split("-")[1].strip()
        if label == "134C":
            if not pd.isna(row[f"Annotated {reporter}"]):  # pyright: ignore[reportGeneralTypeIssues]
                reporter_signal = float(row[f"Annotated {reporter}"])
            if not pd.isna(row["RESGUI_134C Noise"]):  # pyright: ignore[reportGeneralTypeIssues]
                reporter_noise = float(row["RESGUI_134C Noise"])
        else:
            if not pd.isna(row[f"Annotated {reporter}"]):  # pyright: ignore[reportGeneralTypeIssues]
                reporter_signal = float(row[f"Annotated {reporter}"])
            if not pd.isna(row[f"RESGUI_{label} Noise"]):  # pyright: ignore[reportGeneralTypeIssues]
                reporter_noise = float(row[f"RESGUI_{label} Noise"])
        total_signal += reporter_signal
        total_noise += reporter_noise
    return total_signal / total_noise


def __annotate_result_conditions(
    psms: pd.DataFrame, conditions: List[Dict[str, Any]]
) -> pd.DataFrame:
    has_resolution = "RESGUI_Resolution" in psms.columns.tolist()
    if not has_resolution:
        return psms
    for condition in conditions:
        psms[f"Condition_SN_{condition['name']}"] = psms.apply(
            lambda row: __get_sn_for_condition(row, condition["reporters"]), axis=1
        )
    return psms


def __annotate_chimerys_protein_df(
    protein_table: pd.DataFrame,
    psm_table: pd.DataFrame,
    min_chimerys_coefficient: float,
    min_avg_reporter_sn: float,
    min_reporter_res: float,
    min_purity: float,
    conditions: List[Dict[str, Any]],
) -> pd.DataFrame:
    has_resolution = "RESGUI_Resolution" in psm_table.columns.tolist()
    psms_by_proteins = dict()
    for i, psm in tqdm(
        psm_table.iterrows(), total=psm_table.shape[0], desc="Filtering PSMs..."
    ):
        chimerys_coefficient = float(psm["Normalized CHIMERYS Coefficient"])
        avg_reporter_sn = float(psm["Average Reporter S/N"])
        purity = float(psm["Co-Isolation Purity"])
        proteins = [
            protein.strip() for protein in str(psm["Protein Accessions"]).split(";")
        ]
        protein = proteins[0]
        # remove ambiguous PSMs / shared peptides
        if len(proteins) != 1:
            continue
        # remove PSMs with Chimerys Coefficient < threshold
        if (
            pd.isna(chimerys_coefficient)
            or chimerys_coefficient < min_chimerys_coefficient
        ):
            continue
        # remove PSMs with too low average reporter S/N
        if pd.isna(avg_reporter_sn) or avg_reporter_sn < min_avg_reporter_sn:
            continue
        # remove PSMs below purity threshold
        if pd.isna(purity) or purity < min_purity:
            continue
        if protein in psms_by_proteins:
            psms_by_proteins[protein].append(psm)
        else:
            psms_by_proteins[protein] = [psm]
    channels = {key: [] for key in TMT.keys()}
    for i, protein in tqdm(
        protein_table.iterrows(),
        total=protein_table.shape[0],
        desc="Annotating protein abundances...",
    ):
        accession = [x.strip() for x in str(protein["Accession"]).split(";")]
        if len(accession) != 1:
            raise RuntimeError(
                f"Found more then one accession in column 'Accession' in protein table row {i}!"
            )
        accession = accession[0]
        psms_for_accession = list()
        if accession in psms_by_proteins:
            psms_for_accession = psms_by_proteins[accession]
        # else:
        #     print(
        #         f"Info: No PSMs for accession {accession} found due to filter criteria!"
        #     )
        tmt_quants = {key: 0.0 for key in TMT.keys()}
        for psm in psms_for_accession:
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


def __annotate_chimerys_protein_table(
    protein_table: str,
    psm_table: pd.DataFrame,
    settings: Dict[str, Any],
) -> pd.DataFrame:
    protein_df = pd.read_csv(protein_table, sep="\t", low_memory=False)
    min_chimerys_coefficient = float(settings["min_chimerys_coefficient"])
    min_avg_reporter_sn = float(settings["min_avg_reporter_sn"])
    min_reporter_res = float(settings["min_reporter_res"])
    min_purity = float(settings["min_purity"])
    return __annotate_chimerys_protein_df(
        protein_df,
        psm_table,
        min_chimerys_coefficient=min_chimerys_coefficient,
        min_avg_reporter_sn=min_avg_reporter_sn,
        min_reporter_res=min_reporter_res,
        min_purity=min_purity,
        conditions=settings["conditions"],
    )


def __get_bool_from_value(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    elif isinstance(value, int):
        if value in [0, 1]:
            return bool(value)
        else:
            raise ValueError(f"Cannot parse bool value from the given input {value}.")
    elif isinstance(value, str):
        return "t" in value.lower()
    else:
        raise ValueError(f"Cannot parse bool value from the given input {value}.")
    return False


def __get_uncharged_mass_from_exp_mass(mz: float, charge: int) -> float:
    return mz * charge - PROTON * charge


def __get_key(mass: float) -> int:
    return int(round(mass * 10000))


def __check_valid_reporter(reporter: str) -> bool:
    if reporter not in TMT:
        raise ValueError(f"Found invalid reporter label {reporter} in CONDITIONS!")
    return True


def __get_conditions(toml_conditions: Dict[str, Any]) -> List[Dict[str, Any]]:
    conditions = list()
    for cond in toml_conditions["sn_thresholds"]:
        if cond not in toml_conditions:
            raise KeyError(f"Did not find condition {cond} in CONDITIONS!")
        condition = {
            "name": cond,
            "reporters": [
                reporter
                for reporter in toml_conditions[cond]
                if __check_valid_reporter(reporter)
            ],
            "sn": float(toml_conditions["sn_thresholds"][cond]),
        }
        conditions.append(condition)
    return conditions


def __read_settings(toml: str) -> Dict[str, Any]:
    parsed_toml = None
    with open(toml, "rb") as f:
        parsed_toml = tomllib.load(f)
        f.close()
    if parsed_toml is None:
        raise RuntimeError("Could not read config file. Is it in valid TOML format?")
    return {
        "window_size": parsed_toml["METHOD"]["window_size"],
        "window_start": parsed_toml["METHOD"]["window_start"],
        "window_end": parsed_toml["METHOD"]["window_end"],
        "window_overlap": parsed_toml["METHOD"]["window_overlap"],
        "mz_tolerance": parsed_toml["MATCHING"]["mz_tolerance"],
        "rt_tolerance": parsed_toml["MATCHING"]["rt_tolerance"],
        "rt_window": parsed_toml["MATCHING"]["ms1_rt_window"],
        "threshold": parsed_toml["FILTERING"]["total_intensity_threshold"],
        "noise": parsed_toml["FILTERING"]["noise_threshold"],
        "deisotope": parsed_toml["ISOTOPES"]["consider_precursor_isotopes"],
        "isotope_tolerance": parsed_toml["ISOTOPES"]["isotope_tolerance"],
        "max_charge": parsed_toml["ISOTOPES"]["max_charge"],
        "subtract_noise": parsed_toml["QUANTIFICATION"]["subtract_noise"],
        "quantification_method": parsed_toml["QUANTIFICATION"]["quantification_method"],
        "q_value": parsed_toml["PROTEIN"]["q_value"],
        "min_chimerys_coefficient": parsed_toml["PROTEIN"]["min_chimerys_coefficient"],
        "min_avg_reporter_sn": parsed_toml["PROTEIN"]["min_avg_reporter_sn"],
        "min_reporter_res": parsed_toml["PROTEIN"]["min_reporter_res"],
        "min_purity": parsed_toml["PROTEIN"]["min_purity"],
        "keep_pg": parsed_toml["PROTEIN"]["keep_ambiguous_protein_groups"],
        "conditions": __get_conditions(parsed_toml["CONDITIONS"]),
    }


def __get_consensusXML_df(spectrum_filename: str) -> pd.DataFrame:
    in_name = spectrum_filename
    out_name = f"{spectrum_filename}.consensusXML"
    # see https://openms.de/documentation/html/TOPP_IsobaricAnalyzer.html
    subprocess.call(
        [
            "IsobaricAnalyzer",
            "-type",
            "tmt18plex",
            "-in",
            in_name,
            "-out",
            out_name,
            "-ini",
            "tmt18plex_default.ini",
        ]
    )
    # see https://pyopenms.readthedocs.io/en/latest/user_guide/other_ms_data_formats.html#quantiative-data-featurexml-consensusxml
    consensus_features = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(out_name, consensus_features)
    # see https://pyopenms.readthedocs.io/en/latest/user_guide/export_pandas_dataframe.html#consensusmap
    return consensus_features.get_df()


def __get_consensusXML_map(
    consensusXML_df: pd.DataFrame,
) -> Dict[int, Dict[int, pd.Series]]:
    consensusXML_map = dict()
    for i, row in tqdm(
        consensusXML_df.iterrows(),
        total=consensusXML_df.shape[0],
        desc="Reading consensusXML...",
    ):
        mz = __get_key(float(row["mz"]))
        rt = __get_key(float(row["RT"]))
        if mz not in consensusXML_map:
            consensusXML_map[mz] = {rt: row}
        else:
            if rt not in consensusXML_map[mz]:
                consensusXML_map[mz][rt] = row
            else:
                raise RuntimeError(
                    f"Found duplicate entry with m/z {row['mz']} and rt {row['RT']}!"
                )
    return consensusXML_map


def __get_tmt_intensities_oms(
    spectrum: Dict[str, Any], consensusXML_map: Dict[int, Dict[int, pd.Series]]
) -> Dict[str, float]:
    mz_tol = 0.01  # Dalton, m/z
    rt_tol = 1.0  # seconds
    mz = __get_key(spectrum["precursor"])
    rt = __get_key(spectrum["rt"])
    mz_tol = __get_key(mz_tol)
    rt_tol = __get_key(rt_tol)
    row = None
    for i in range(mz_tol):
        if (mz - i) in consensusXML_map:
            for j in range(rt_tol):
                if (rt - j) in consensusXML_map[mz - i]:
                    row = consensusXML_map[mz - i][rt - j]
                    break
                if (rt + j) in consensusXML_map[mz - i]:
                    row = consensusXML_map[mz - i][rt + j]
                    break
            if row is not None:
                break
        if (mz + i) in consensusXML_map:
            for j in range(rt_tol):
                if (rt - j) in consensusXML_map[mz + i]:
                    row = consensusXML_map[mz + i][rt - j]
                    break
                if (rt + j) in consensusXML_map[mz + i]:
                    row = consensusXML_map[mz + i][rt + j]
                    break
            if row is not None:
                break
    tmt_quants = {key: 0.0 for key in TMT.keys()}
    if row is None:
        return tmt_quants
    for k, v in TMT_OMS.items():
        tmt_quants[k] += row[v]
    return tmt_quants


def __get_tmt_intensities(spectrum: Dict[str, Any]) -> Dict[str, float]:
    # this does not do any kind of isotope corrections, for that purpose use
    # __get_tmt_intensities_oms
    # instead!
    tmt_quants = {key: 0.0 for key in TMT.keys()}
    for reporter_ion_name, reporter_ion_mass in TMT.items():
        for i, mz in enumerate(spectrum["mz_array"]):
            if __within_tolerance(mz, reporter_ion_mass, TMT_TOLERANCE):
                tmt_quants[reporter_ion_name] += spectrum["intensity_array"][i]
                break
    return tmt_quants


def __get_resolution_gui_map(filename: str) -> Dict[str, Dict[int, pd.Series]]:
    resolution_gui_map = dict()
    df = pd.read_csv(filename, low_memory=False)
    for i, row in tqdm(
        df.iterrows(), total=df.shape[0], desc="Reading Resolution GUI output..."
    ):
        _head, tail = os.path.split(str(row["Raw File"]).strip())
        spectrum_filename = tail.strip()
        spectrum_scannr = int(row["MS2/MS3 Scan"])
        if spectrum_filename in resolution_gui_map:
            if spectrum_scannr in resolution_gui_map[spectrum_filename]:
                warnings.warn(
                    RuntimeWarning(
                        f"Found duplicate MS2 scan with scan number {spectrum_scannr} in Resolution GUI output! Using first MS2 scan!"
                    )
                )
            else:
                resolution_gui_map[spectrum_filename][spectrum_scannr] = row
        else:
            resolution_gui_map[spectrum_filename] = {spectrum_scannr: row}
    return resolution_gui_map


def __get_resolution_gui_values(
    spectrum: Dict[str, Any],
    resolution_gui_map: Dict[str, Dict[int, pd.Series]],
    spectrum_filename: str,
) -> Dict[str, float]:
    _head, tail = os.path.split(spectrum_filename)
    sf = tail[:-5] + ".raw" if tail[-5:].lower() == ".mzml" else tail
    sn = spectrum["scan_nr"]
    row = resolution_gui_map[sf][sn]
    data = dict()
    for col in RESOLUTION_GUI_COLS:
        data[f"RESGUI_{col}"] = row[col]
    return data


def __get_tmt_intensities_resgui(
    spectrum: Dict[str, Any],
    resolution_gui_map: Dict[str, Dict[int, pd.Series]],
    spectrum_filename: str,
) -> Dict[str, float]:
    resolution_gui_values = __get_resolution_gui_values(
        spectrum, resolution_gui_map, spectrum_filename
    )
    tmt_quants = {key: 0.0 for key in TMT.keys()}
    for reporter_ion_name in TMT.keys():
        channel = reporter_ion_name.split("-")[1].strip()
        descriptor = "Intesntiy" if channel == "134C" else "Intensity"
        resgui_col = f"RESGUI_{channel} {descriptor}"
        tmt_quants[reporter_ion_name] += resolution_gui_values[resgui_col]
    return tmt_quants


def __subtract_noise(
    tmt_quants: Dict[str, float],
    spectrum: Dict[str, Any],
    resolution_gui_map: Dict[str, Dict[int, pd.Series]],
    spectrum_filename: str,
) -> Dict[str, float]:
    resolution_gui_values = __get_resolution_gui_values(
        spectrum, resolution_gui_map, spectrum_filename
    )
    new_tmt_quants = dict()
    for reporter_ion_name in TMT.keys():
        channel = reporter_ion_name.split("-")[1].strip()
        resgui_col = f"RESGUI_{channel} Noise"
        new_tmt_quants[reporter_ion_name] = (
            tmt_quants[reporter_ion_name] - resolution_gui_values[resgui_col]
        )
    return new_tmt_quants


def __parse_scan_nr_from_id(id: str) -> int:
    return int(id.split("scan=")[1].split()[0])


def __read_spectra_by_scannumber(
    filename: str,
) -> Dict[str, Dict[int, Dict[str, Any]]]:
    spectra_ms1 = dict()
    spectra_ms2 = dict()
    total = 0
    total_ms1 = 0
    total_ms2 = 0
    with mzml.read(filename) as reader:
        for spectrum in reader:
            scan = __parse_scan_nr_from_id(spectrum["id"])
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
                        s = dict()
                        s["scan_nr"] = scan
                        s["precursor"] = float(ion["selected ion m/z"])
                        s["rt"] = rt_in_sec
                        s["mz_array"] = spectrum["m/z array"]
                        s["intensity_array"] = spectrum["intensity array"]
                        if scan in spectra_ms2:
                            raise RuntimeError()
                        else:
                            spectra_ms2[scan] = s
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


# checks if two values (value and ref) are equal with a certain tolerance (tol)
def __within_tolerance(value: float, ref: float, tol: float) -> bool:
    return (value > ref - tol) and (value < ref + tol)


# checks if the specified precursor m/z can be found in the given MS1 spectrum
# considering tolerance
def __check_mz_in_ms1(
    precursor_mz: float, spectrum: Dict[str, Any], mz_tol: float
) -> bool:
    for mz in spectrum["mz_array"]:
        if __within_tolerance(mz, precursor_mz, mz_tol):
            return True
    return False


def __calculate_purity(
    spectrum: Dict[str, Any],  # MS1 spectrum
    matching_window: Tuple[float, float],
    precursor: float,
    precursor_index: int,
    precursor_intensity: float,
    max_charge: int,
    noise_threshold: float,
    do_deisotope: bool,
    isotope_tol: float,
) -> float:
    # get highest intensity peak in window
    most_intense_peak = 0.0
    peaks_in_window_mz = list()
    peaks_in_window_in = list()
    for i in range(len(spectrum["mz_array"])):
        if (
            spectrum["mz_array"][i] > matching_window[0]
            and spectrum["mz_array"][i] <= matching_window[1]
        ):
            if spectrum["intensity_array"][i] > most_intense_peak:
                most_intense_peak = spectrum["intensity_array"][i]
            peaks_in_window_mz.append(spectrum["mz_array"][i])
            peaks_in_window_in.append(spectrum["intensity_array"][i])
    if do_deisotope:
        # look for precursor isotope peaks
        charge = None
        ## determine charge
        for i in range(1, 11):  # +- 10 peaks are considered for charge determination
            if precursor_index + i < len(spectrum["mz_array"]):
                curr_peak_mz = spectrum["mz_array"][precursor_index + i]
                for c in range(1, max_charge + 1):
                    if __within_tolerance(
                        curr_peak_mz, precursor + ISOTOPE / c, isotope_tol
                    ):
                        charge = c
                        break
            if charge is None:
                if precursor_index - i >= 0:
                    curr_peak_mz = spectrum["mz_array"][precursor_index - i]
                    for c in range(1, max_charge + 1):
                        if __within_tolerance(
                            curr_peak_mz, precursor - ISOTOPE / c, isotope_tol
                        ):
                            charge = c
                            break
            if charge is not None:
                break
        # get precursor isotope peaks
        if charge is not None:
            for i in range(1, 6):  # +- 5 isotopes are considered for every precursor
                precursor_isotope_plus_i = precursor + i * (ISOTOPE / charge)
                if precursor_isotope_plus_i <= matching_window[1]:
                    for j in range(len(peaks_in_window_mz)):
                        if peaks_in_window_in[j] / most_intense_peak >= noise_threshold:
                            if __within_tolerance(
                                peaks_in_window_mz[j],
                                precursor_isotope_plus_i,
                                isotope_tol,
                            ):
                                precursor_intensity += peaks_in_window_in[j]
                                break
                precursor_isotope_minus_i = precursor - i * (ISOTOPE / charge)
                if precursor_isotope_minus_i > matching_window[0]:
                    for j in range(len(peaks_in_window_mz)):
                        if peaks_in_window_in[j] / most_intense_peak >= noise_threshold:
                            if __within_tolerance(
                                peaks_in_window_mz[j],
                                precursor_isotope_minus_i,
                                isotope_tol,
                            ):
                                precursor_intensity += peaks_in_window_in[j]
                                break
    # if precursor is noisy, return 0.0
    if precursor_intensity / most_intense_peak < noise_threshold:
        return 0.0
    # calculate total intensity in window
    total_intensity_in_window = 0.0
    for i in range(len(spectrum["mz_array"])):
        if (
            spectrum["mz_array"][i] > matching_window[0]
            and spectrum["mz_array"][i] <= matching_window[1]
        ):
            if spectrum["intensity_array"][i] / most_intense_peak > noise_threshold:
                total_intensity_in_window += spectrum["intensity_array"][i]
    # return intensity ratio (purity)
    return precursor_intensity / total_intensity_in_window


def __get_optimal_window(
    windows: List[Tuple[float, float]], spectrum: Dict[str, Any]
) -> List[Tuple[float, float]]:
    deviations = list()
    for window in windows:
        center = window[0] + (window[1] - window[0]) / 2.0
        deviation = abs(center - spectrum["precursor"])
        deviations.append(deviation)
    return [windows[deviations.index(min(deviations))]]


# calculates the purity for a given precursor
def __calculate_precursor_intensity_ms1(
    precursor_mz: float,
    spectrum: Dict[str, Any],  # MS1 spectrum
    spectrum_ms2: Dict[str, Any],  # MS2 spectrum
    mz_tol: float,
    do_deisotope: bool,
    isotope_tol: float,
    max_charge: int,
    noise_threshold: float,
    windows: List[Tuple[float, float]],
) -> float | None:
    # parameter windows should define all DIA windows including their start and
    # end points (in m/z)
    precursor = None
    precursor_index = None
    precursor_intensity = None
    # first look for the precursor in the MS1 spectrum
    for i in range(len(spectrum["mz_array"])):
        if __within_tolerance(spectrum["mz_array"][i], precursor_mz, mz_tol):
            # if no precursor yet found, set precursor
            if precursor is None:
                precursor = spectrum["mz_array"][i]
                precursor_index = i
                precursor_intensity = spectrum["intensity_array"][i]
            # else we need to deal with the fact that there is multiple potential
            # precursors
            else:
                # todo clarify behaviour
                if STRATEGY == 1:
                    # use precursor with highest intensity
                    if spectrum["intensity_array"][i] > precursor_intensity:
                        precursor = spectrum["mz_array"][i]
                        precursor_index = i
                        precursor_intensity = spectrum["intensity_array"][i]
                elif STRATEGY == 2:
                    # do not use identification
                    return None
                elif STRATEGY == 3:
                    # use closest precursor
                    raise NotImplementedError()
                else:
                    raise RuntimeError(
                        f"Found ambiguous precursors in MS1 spectrum for precursor m/z {precursor_mz} using MS1 spectrum at retention time {spectrum['rt']}."
                    )
    # if no precursor is found an error is raised
    if precursor is None or precursor_index is None or precursor_intensity is None:
        warnings.warn(
            RuntimeWarning(
                f"Could not find precursor in MS1 spectrum with scan number {spectrum['scan_nr']}."
            )
        )
        return None
    # find the corresponding m/z window(s) that the precursor is in
    matching_windows = list()
    for window in windows:
        if precursor > window[0] and precursor <= window[1]:
            matching_windows.append(window)
    # if no matching window is found an error is raised
    if len(matching_windows) == 0:
        raise RuntimeError(
            f"Could not find matching window for precursor m/z {precursor_mz}!"
        )
    # select optimal window?
    select_optimal_window = True
    if len(matching_windows) > 1 and select_optimal_window:
        matching_windows = __get_optimal_window(matching_windows, spectrum_ms2)
    # if there are multiple matching windows, we calculate the purity of every window
    # and return the minimum purity
    purities = list()
    for matching_window in matching_windows:
        purities.append(
            __calculate_purity(
                spectrum,
                matching_window,
                precursor,
                precursor_index,
                precursor_intensity,
                max_charge,
                noise_threshold,
                do_deisotope,
                isotope_tol,
            )
        )
    return min(purities)


def __get_ms2_spectrum_by_scannumber(
    scan_nr: int,
    precursor_mz: float,
    mz_tol: float,
    retention_time_ms1_window: float,
    do_deisotope: bool,
    isotope_tol: float,
    max_charge: int,
    noise_threshold: float,
    spectra: Dict[str, Any],
    windows: List[Tuple[float, float]],
) -> Dict[str, Any]:
    spectrum = spectra["ms2"][scan_nr]
    retention_time_key = __get_key(spectrum["rt"])
    # look for closest MS1 spectrum that has precursor (with m/z tolerance mz_tol)
    # within ms1 rt window
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
                f"Could not find a suitable MS1 spectrum for precursor m/z {precursor_mz} and retention time {spectrum['rt']}. MS2 scan: {scan_nr}"
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


# given window start and end points (in m/z) and a window size, calculate all
# m/z windows
def __get_windows(
    window_start: float,
    window_end: float,
    window_size: float,
    window_overlap: float,
    window_file: Optional[str],
) -> List[Tuple[float, float]]:
    windows = list()
    if window_file is None:
        current_window_start = window_start
        while current_window_start < window_end:
            current_window_end = current_window_start + window_size
            if current_window_end > window_end:
                windows.append((current_window_start, window_end))
                break
            windows.append((current_window_start, current_window_end))
            current_window_start = current_window_start + window_size - window_overlap
    else:
        df = pd.read_csv(window_file)
        for i, row in df.iterrows():
            w1 = str(row["m/z range"]).split("-")[0]
            w2 = str(row["m/z range"]).split("-")[1]
            windows.append((float(w1), float(w2)))
    return windows


# annotates the Chimerys result with purity and TMT quantities
# currently based on the PSM table
def __annotate_chimerys_result(
    filename: str,
    spectrum_filename: str,
    spectra: Dict[str, Any],
    settings: Dict[str, Any],
    consensusXML_map: Optional[Dict[int, Dict[int, pd.Series]]] = None,
    resolution_gui_map: Optional[Dict[str, Dict[int, pd.Series]]] = None,
    window_file: Optional[str] = None,
) -> pd.DataFrame:
    # spectra should be given by __read_spectra_by_scannumber
    # settings should be given by __read_settings
    df = pd.read_csv(filename, sep="\t", low_memory=False)
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
        df.iterrows(), total=df.shape[0], desc="Annotating Chimerys result..."
    ):
        # get scan number
        scan_nr = int(row["First Scan"])
        # get Chimerys precursor m/z from identification
        prec_mz = float(row["m/z [Da]"])
        # settings defined m/z tolerance
        mz_tol = float(settings["mz_tolerance"])
        # settings definied retention time window for MS1 spectra
        rt_window = float(settings["rt_window"])
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
        spectrum_purity = __get_ms2_spectrum_by_scannumber(
            scan_nr,
            prec_mz,
            mz_tol,
            rt_window,
            do_deisotope,
            isotope_tolerance,
            max_charge,
            noise_threshold,
            spectra,
            windows,
        )
        spectrum = spectrum_purity["spectrum"]
        purity = spectrum_purity["purity"]
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
    # update Chimerys result
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


@Gooey(
    encoding="utf-8",
    program_name=f"TMT Chimerys DIA {__version}",
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
        prog="tmt_chimerys.py",
        description="Calculates co-isolation purity for Chimerys DIA TMT PSMs and optionally quantifies them.",
        epilog="(c) Research Institute of Molecular Pathology, 2025",
    )
    req = parser.add_argument_group("Required", "Required Arguments.")
    req.add_argument(
        "-i",
        "--chimerys",
        dest="chimerys",
        required=True,
        help="Path/name of the Chimerys PSM result file in tab-separated .txt format.",
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
        required=True,
        help="Path/name of the config file.",
        type=str,
        widget="FileChooser",
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
        widget="FileChooser",
    )
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
        dest="window_file",
        default=None,
        help="Window file, overrides config file!",
        type=str,
        widget="FileChooser",
    )
    args = parser.parse_args(argv)
    settings = __read_settings(args.config)
    print("Read settings:")
    print(settings)
    if args.window_file is not None:
        print(f"Using windows from given windows file: {args.window_file}")
    args_spectra = __convert(args.spectra)
    spectra = __read_spectra_by_scannumber(args_spectra)
    quantification_method = int(settings["quantification_method"])
    consensusXML_map = None
    if quantification_method != 1 and quantification_method != 3:
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
        window_file=args.window_file,
    )
    df.to_csv(
        args.chimerys.split(".txt")[0] + "_purity_tmt_quant.txt",
        sep="\t",
        index=False,
    )
    df = __annotate_result_conditions(df, settings["conditions"])
    df.to_csv(
        args.chimerys.split(".txt")[0] + "_purity_tmt_quant_conditions.txt",
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
