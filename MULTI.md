# Running the Scripts for Multiple Files

If you want to run the `tmt_chimerys.py` or `tmt_chimerys_dda.py` script for multiple inputs, please
install the requirements and follow the steps below.

## Requirements

- Please install [OpenMS](https://openms.readthedocs.io/en/latest/about/installation.html).
- Please install either [python 3.12 or higher](https://www.python.org/downloads/)
  or [uv](https://docs.astral.sh/uv/).

## Steps

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
> Please note that the same config file is used for all the files in the folder. The DIA script automatically
> parses the window size from the MS file name and overrides the window size in the config file. You might have
> to adjust this parsing procedure if your filenames do not follow the same naming pattern.
