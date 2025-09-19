# Running the Scripts for Multiple Files

If you want to run the `tmt_chimerys.py` or `tmt_chimerys_dda.py` script for multiple inputs, please
follow the following steps:

- Put all your `.raw` files and identification files (e.g. PSMs and optionally proteins) in the same folder.
  - For example, let's call this folder `tmt_files`.
- PSM files should end in the suffix `_PSMs.txt` and be in tab-separated format.
- Protein files should end in the suffix `_Proteins.txt` and be in tab-separated format.
- For example, like this:
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1.raw`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1_PSMs.txt`
  - `20250519_Astral1_Evo_TH070_TT_THIDmulti003_pool_DIA_mz0c5_3ng_1_Proteins.txt`
- It's important that files that belong together are named like this schema.
- It's not required to convert `.raw` files to `.mzML` format but it substantially speeds up the process.
  Here's how to do it:
  - Download ThermoRawFileParser from [here](https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v1.4.5).
  - Convert your RAW files with:
    ```bash
    ThermoRawFileParser.exe -d path/to/tmt_files
    ```
- Running the multi-file scripts also requires that the `resolution.csv` from the Resolution GUI tool is in the same.
- Additionally, your `config.toml` should also be in the same folder. Adapt the configuration to your needs.
- Lastly, please put the `tmt18plex_default.ini`
