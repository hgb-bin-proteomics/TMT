# Commandline Interface

You can run the following scripts from the commandline using [python](https://www.python.org/downloads/)
or [uv](https://docs.astral.sh/uv/).

You can find all scripts in the `/scripts` folder.

> [!IMPORTANT]
>
> Please make sure that `tmt_chimerys.py`, `tmt_chimerys_dda.py`, `tmt_diann.py`, `tmt_spectronaut.py`,
> and `tmt18plex_default.ini` are in the same directory when running them as python scripts!

## Chimerys DIA

- Export Chimerys PSMs from Proteome Discoverer in tab-separated `.txt` format.
- \[Optionally\] Export Chimerys Proteins from Proteome Discoverer in tab-seperated `.txt` format.
- Set you desired parameters in `config.toml`.
- The scripts support both `.raw` files and `.mzML` files as input, `.raw` files will be automatically
  converted to `.mzML` when the scripts are run.
- The following steps are optional if you want to convert your `.raw` files manually:
  - Download ThermoRawFileParser from [here](https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v1.4.5).
  - Convert your RAW file with:
    ```bash
    ThermoRawFileParser.exe -i RAW_FILE_NAME.raw
    ```
- Install [OpenMS](https://openms.readthedocs.io/en/latest/about/installation.html).
- **Option A (recommended): Running via [uv](https://docs.astral.sh/uv/).**
  - [Install uv](https://docs.astral.sh/uv/getting-started/installation/) if it's not already installed on your system, e.g.:
    ```bash
    pip install uv
    ```
  - Run the script with:
    ```bash
    uv run tmt_chimerys.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMs.txt -c config.toml
    ```
  - _or_ if you also have proteins with:
    ```bash
    uv run tmt_chimerys.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMs.txt -c config.toml -p PROTEOME_DISCOVERER_Proteins.txt
    ```
  - To display all possible parameters run:
    ```bash
    uv run tmt_chimerys.py --help
    ```
  - Alternatively you can also run the script with a graphical user interface using:
    ```bash
    uv run tmt_chimerys_gui.py
    ```
- **Option B: Running via native python.**
  - Install python 3.12 or greater from [here](https://www.python.org/downloads/).
  - Install requirements with:
    ```bash
    pip install -r requirements.txt
    ```
  - Run the script with:
    ```bash
    python tmt_chimerys.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMs.txt -c config.toml
    ```
  - _or_ if you also have proteins with:
    ```bash
    python tmt_chimerys.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMs.txt -c config.toml -p PROTEOME_DISCOVERER_Proteins.txt
    ```
  - To display all possible parameters run:
    ```bash
    python tmt_chimerys.py --help
    ```
  - Alternatively you can also run the script with a graphical user interface using:
    ```bash
    python tmt_chimerys_gui.py
    ```
- The result will be new files with name extension `_purity_tmt_quant` that are written out,
  containing purity and quantification values.

## Chimerys DDA

- Export Chimerys PSMs from Proteome Discoverer in tab-separated `.txt` format.
- \[Optionally\] Export Chimerys Proteins from Proteome Discoverer in tab-seperated `.txt` format.
- Set you desired parameters in `config.toml`.
- The scripts support both `.raw` files and `.mzML` files as input, `.raw` files will be automatically
  converted to `.mzML` when the scripts are run.
- The following steps are optional if you want to convert your `.raw` files manually:
  - Download ThermoRawFileParser from [here](https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v1.4.5)
  - Convert your RAW file with:
    ```bash
    ThermoRawFileParser.exe -i RAW_FILE_NAME.raw
    ```
- Install [OpenMS](https://openms.readthedocs.io/en/latest/about/installation.html).
- **Option A (recommended): Running via [uv](https://docs.astral.sh/uv/).**
  - [Install uv](https://docs.astral.sh/uv/getting-started/installation/) if it's not already installed on your system, e.g.:
    ```bash
    pip install uv
    ```
  - Run the script with:
    ```bash
    uv run tmt_chimerys_dda.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMs.txt -c config.toml
    ```
  - _or_ if you also have proteins with:
    ```bash
    uv run tmt_chimerys_dda.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMs.txt -c config.toml -p PROTEOME_DISCOVERER_Proteins.txt
    ```
  - To display all possible parameters run:
    ```bash
    uv run tmt_chimerys_dda.py --help
    ```
  - Alternatively you can also run the script with a graphical user interface using:
    ```bash
    uv run tmt_chimerys_dda_gui.py
    ```
- **Option B: Running via native python.**
  - Install python 3.12 or greater from [here](https://www.python.org/downloads/).
  - Install requirements with:
    ```bash
    pip install -r requirements.txt
    ```
  - Run the script with:
    ```bash
    python tmt_chimerys_dda.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMs.txt -c config.toml
    ```
  - _or_ if you also have proteins with:
    ```bash
    python tmt_chimerys_dda.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMs.txt -c config.toml -p PROTEOME_DISCOVERER_Proteins.txt
    ```
  - To display all possible parameters run:
    ```bash
    python tmt_chimerys_dda.py --help
    ```
  - Alternatively you can also run the script with a graphical user interface using:
    ```bash
    python tmt_chimerys_dda_gui.py
    ```
- The result will be new files with name extension `_purity_tmt_quant` that are written out,
  containing purity and quantification values.

## Spectronaut

- Export matched precursors/the main report from Spectronaut in comma-separated `.csv` format.
- Set you desired parameters in `config.toml`.
- The scripts support both `.raw` files and `.mzML` files as input, `.raw` files will be automatically
  converted to `.mzML` when the scripts are run.
- The following steps are optional if you want to convert your `.raw` files manually:
  - Download ThermoRawFileParser from [here](https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v1.4.5)
  - Convert your RAW file with:
    ```bash
    ThermoRawFileParser.exe -i RAW_FILE_NAME.raw
    ```
- Install [OpenMS](https://openms.readthedocs.io/en/latest/about/installation.html).
- **Option A (recommended): Running via [uv](https://docs.astral.sh/uv/).**
  - [Install uv](https://docs.astral.sh/uv/getting-started/installation/) if it's not already installed on your system, e.g.:
    ```bash
    pip install uv
    ```
  - Run the script with:
    ```bash
    uv run tmt_spectronaut.py -s SPECTRA.mzML -i report.csv -c config.toml
    ```
  - To display all possible parameters run:
    ```bash
    uv run tmt_spectronaut.py --help
    ```
  - Alternatively you can also run the script with a graphical user interface using:
    ```bash
    uv run tmt_spectronaut_gui.py
    ```
- **Option B: Running via native python.**
  - Install python 3.12 or greater from [here](https://www.python.org/downloads/).
  - Install requirements with:
    ```bash
    pip install -r requirements.txt
    ```
  - Run the script with:
    ```bash
    python tmt_spectronaut.py -s SPECTRA.mzML -i report.csv -c config.toml
    ```
  - To display all possible parameters run:
    ```bash
    python tmt_spectronaut.py --help
    ```
  - Alternatively you can also run the script with a graphical user interface using:
    ```bash
    python tmt_spectronaut_gui.py
    ```
- The result will be new files with name extension `_purity_tmt_quant` that are written out,
  containing purity and quantification values.

## DIA-NN

- Use the `report.parquet` that you get from DIA-NN.
- Set you desired parameters in `config.toml`.
- The scripts support both `.raw` files and `.mzML` files as input, `.raw` files will be automatically
  converted to `.mzML` when the scripts are run.
- The following steps are optional if you want to convert your `.raw` files manually:
  - Download ThermoRawFileParser from [here](https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v1.4.5)
  - Convert your RAW file with:
    ```bash
    ThermoRawFileParser.exe -i RAW_FILE_NAME.raw
    ```
- Install [OpenMS](https://openms.readthedocs.io/en/latest/about/installation.html).
- **Option A (recommended): Running via [uv](https://docs.astral.sh/uv/).**
  - [Install uv](https://docs.astral.sh/uv/getting-started/installation/) if it's not already installed on your system, e.g.:
    ```bash
    pip install uv
    ```
  - Run the script with:
    ```bash
    uv run tmt_diann.py -s SPECTRA.mzML -i report.parquet -c config.toml
    ```
  - To display all possible parameters run:
    ```bash
    uv run tmt_diann.py --help
    ```
  - Alternatively you can also run the script with a graphical user interface using:
    ```bash
    uv run tmt_diann_gui.py
    ```
- **Option B: Running via native python.**
  - Install python 3.12 or greater from [here](https://www.python.org/downloads/).
  - Install requirements with:
    ```bash
    pip install -r requirements.txt
    ```
  - Run the script with:
    ```bash
    python tmt_diann.py -s SPECTRA.mzML -i report.parquet -c config.toml
    ```
  - To display all possible parameters run:
    ```bash
    python tmt_diann.py --help
    ```
  - Alternatively you can also run the script with a graphical user interface using:
    ```bash
    python tmt_diann_gui.py
    ```
- The result will be new files with name extension `_purity_tmt_quant` that are written out,
  containing purity and quantification values.
