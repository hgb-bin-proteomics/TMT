# Running the Scripts for Multiple Files

If you want to run the scripts for multiple input files sequentially, please
install the requirements and follow the steps below.

## Requirements

- Please install [OpenMS](https://openms.readthedocs.io/en/latest/about/installation.html).
- Please install either [python 3.12 or higher](https://www.python.org/downloads/)
  or [uv](https://docs.astral.sh/uv/).

## Chimerys: Steps

- Put all your `.raw` files and identification files (e.g. PSMs and optionally proteins) in the same folder.
  - For example, let's call this folder `tmt_files`.
- PSM files should end in the suffix `_PSMs.txt` and be in tab-separated format.
- Protein files should end in the suffix `_Proteins.txt` and be in tab-separated format.
- For example, like this:
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.raw`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1_PSMs.txt`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1_Proteins.txt`
- It's important that files that belong together are named like this schema.
- You will need to convert your `.raw` files to `.mzML` format, here's how to do it:
  - Download ThermoRawFileParser from [here](https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v1.4.5).
  - Convert your RAW files with:
    ```bash
    ThermoRawFileParser.exe -d path/to/tmt_files
    ```
- Running the multi-file scripts also requires that the `resolution.csv` from the
  [Resolution GUI tool](https://doi.org/10.1038/s41467-025-60022-x) is in the same folder.
- Additionally, your `config.toml` should also be in the same folder. Adapt the configuration to your needs.
- Moreover, please put the `tmt18plex_default.ini` also in that folder and adapt if needed
  (see [here](https://openms.de/documentation/html/TOPP_IsobaricAnalyzer.html)).
- If you have a file with DIA mass windows, you should also put it in that folder.
  - You also have to set the filename in the `tmt_chimerys_multi.py` script under `WINDOW_FILE`.
  - You can ignore this for DDA results.
- Lastly, please put the following scripts into the folder:
  - `scripts/tmt_chimerys.py`
  - `scripts/tmt_chimerys_dda.py`
  - `multi/tmt_chimerys_multi.py`
  - `multi/tmt_chimerys_dda_multi.py`
- Your `tmt_files` folder structure should now be something like this:
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.raw`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.mzML`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1_PSMs.txt`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1_Proteins.txt`
  - `...`
  - `resolution.csv`
  - `config.toml`
  - `tmt18plex_default.ini`
  - `Mass List Table.csv` _(optional)_
  - `tmt_chimerys.py`
  - `tmt_chimerys_dda.py`
  - `tmt_chimerys_multi.py`
  - `tmt_chimerys_dda_multi.py`
- Open a terminal in this folder.
- **Option A (recommended): Running via [uv](https://docs.astral.sh/uv/).**
  - [Install uv](https://docs.astral.sh/uv/getting-started/installation/) if it's not already installed on your system, e.g.:
    ```bash
    pip install uv
    ```
  - Run the DIA script with:
    ```bash
    uv run tmt_chimerys_multi.py
    ```
  - _or_ run the DDA script with:
    ```bash
    uv run tmt_chimerys_dda_multi.py
    ```
- **Option B: Running via native python.**
  - Install python 3.12 or greater from [here](https://www.python.org/downloads/).
  - Install requirements from the `requirements.txt` file in the root directory with:
    ```bash
    pip install -r requirements.txt
    ```
  - Run the DIA script with:
    ```bash
    python tmt_chimerys_multi.py
    ```
  - _or_ run the DDA script with:
    ```bash
    python tmt_chimerys_dda_multi.py
    ```

> [!IMPORTANT]
>
> Please note that the same config file is used for all the files in the folder. The scripts automatically
> parse the window size from the MS file name and overrides the window size in the config file. You might have
> to adjust this parsing procedure if your filenames do not follow the same naming pattern.
> If a `WINDOW_FILE` is given, it will always use the windows from the file!

## DIA-NN: Steps

- Put all your `.raw` files and the identification file (e.g. precursors/the main report) in the same folder.
  - For example, let's call this folder `tmt_files`.
- The report file should end in the suffix `.parquet` and be in parquet format.
- For example, like this:
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.raw`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.parquet`
- Or like this:
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.raw`
  - `report.parquet`
- You need to set the filename of the identification file in the `tmt_diann_multi.py` file, e.g. `MAIN_REPORT="report.parquet"`.
- You will need to convert your `.raw` files to `.mzML` format, here's how to do it:
  - Download ThermoRawFileParser from [here](https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v1.4.5).
  - Convert your RAW files with:
    ```bash
    ThermoRawFileParser.exe -d path/to/tmt_files
    ```
- Running the multi-file scripts also requires that the `resolution.csv` from the
  [Resolution GUI tool](https://doi.org/10.1038/s41467-025-60022-x) is in the same folder.
- Additionally, your `config.toml` should also be in the same folder. Adapt the configuration to your needs.
- Moreover, please put the `tmt18plex_default.ini` also in that folder and adapt if needed
  (see [here](https://openms.de/documentation/html/TOPP_IsobaricAnalyzer.html)).
- If you have a file with DIA mass windows, you should also put it in that folder.
  - You also have to set the filename in the `tmt_chimerys_multi.py` script under `WINDOW_FILE`.
  - You can ignore this for DDA results.
- Lastly, please put the following scripts into the folder:
  - `scripts/tmt_chimerys.py`
  - `scripts/tmt_spectronaut.py`
  - `scripts/tmt_diann.py`
  - `multi/tmt_diann_multi.py`
- Your `tmt_files` folder structure should now be something like this:
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.raw`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.mzML`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1_precursor.parquet`
  - `...`
  - `resolution.csv`
  - `config.toml`
  - `tmt18plex_default.ini`
  - `Mass List Table.csv` _(optional)_
  - `tmt_chimerys.py`
  - `tmt_spectronaut.py`
  - `tmt_diann.py`
  - `tmt_diann_multi.py`
- Open a terminal in this folder.
- **Option A (recommended): Running via [uv](https://docs.astral.sh/uv/).**
  - [Install uv](https://docs.astral.sh/uv/getting-started/installation/) if it's not already installed on your system, e.g.:
    ```bash
    pip install uv
    ```
  - Run the DIA-NN script with:
    ```bash
    uv run tmt_diann_multi.py
    ```
- **Option B: Running via native python.**
  - Install python 3.12 or greater from [here](https://www.python.org/downloads/).
  - Install requirements from the `requirements.txt` file in the root directory with:
    ```bash
    pip install -r requirements.txt
    ```
  - Run the DIA-NN script with:
    ```bash
    python tmt_diann_multi.py
    ```

> [!IMPORTANT]
>
> Please note that the same config file is used for all the files in the folder. The scripts automatically
> parse the window size from the MS file name and overrides the window size in the config file. You might have
> to adjust this parsing procedure if your filenames do not follow the same naming pattern.
> If a `WINDOW_FILE` is given, it will always use the windows from the file!

## Spectronaut: Steps

- Put all your `.raw` files and identification file (e.g. precursors/the main report) in the same folder.
  - For example, let's call this folder `tmt_files`.
- The report file should end in the suffix `.csv` and be in comma-separated (`.csv`) format.
- For example, like this:
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.raw`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.csv`
- Or like this:
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.raw`
  - `report.csv`
- You need to set the filename of the identification file in the `tmt_spectronaut_multi.py` file, e.g. `MAIN_REPORT="report.csv"`.
- You will need to convert your `.raw` files to `.mzML` format, here's how to do it:
  - Download ThermoRawFileParser from [here](https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v1.4.5).
  - Convert your RAW files with:
    ```bash
    ThermoRawFileParser.exe -d path/to/tmt_files
    ```
- Running the multi-file scripts also requires that the `resolution.csv` from the
  [Resolution GUI tool](https://doi.org/10.1038/s41467-025-60022-x) is in the same folder.
- Additionally, your `config.toml` should also be in the same folder. Adapt the configuration to your needs.
- Moreover, please put the `tmt18plex_default.ini` also in that folder and adapt if needed
  (see [here](https://openms.de/documentation/html/TOPP_IsobaricAnalyzer.html)).
- If you have a file with DIA mass windows, you should also put it in that folder.
  - You also have to set the filename in the `tmt_chimerys_multi.py` script under `WINDOW_FILE`.
  - You can ignore this for DDA results.
- Lastly, please put the following scripts into the folder:
  - `scripts/tmt_chimerys.py`
  - `scripts/tmt_spectronaut.py`
  - `multi/tmt_spectronaut_multi.py`
- Your `tmt_files` folder structure should now be something like this:
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.raw`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.mzML`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1_precursor.csv`
  - `...`
  - `resolution.csv`
  - `config.toml`
  - `tmt18plex_default.ini`
  - `Mass List Table.csv` _(optional)_
  - `tmt_chimerys.py`
  - `tmt_spectronaut.py`
  - `tmt_spectronaut_multi.py`
- Open a terminal in this folder.
- **Option A (recommended): Running via [uv](https://docs.astral.sh/uv/).**
  - [Install uv](https://docs.astral.sh/uv/getting-started/installation/) if it's not already installed on your system, e.g.:
    ```bash
    pip install uv
    ```
  - Run the Spectronaut script with:
    ```bash
    uv run tmt_spectronaut_multi.py
    ```
- **Option B: Running via native python.**
  - Install python 3.12 or greater from [here](https://www.python.org/downloads/).
  - Install requirements from the `requirements.txt` file in the root directory with:
    ```bash
    pip install -r requirements.txt
    ```
  - Run the Spectronaut script with:
    ```bash
    python tmt_spectronaut_multi.py
    ```

> [!IMPORTANT]
>
> Please note that the same config file is used for all the files in the folder. The scripts automatically
> parse the window size from the MS file name and overrides the window size in the config file. You might have
> to adjust this parsing procedure if your filenames do not follow the same naming pattern.
> If a `WINDOW_FILE` is given, it will always use the windows from the file!
