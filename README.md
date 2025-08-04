# TMT

## Usage

> [!IMPORTANT]
> Please make sure that `tmt_chimerys.py`, `tmt_chimerys_dda.py`, `tmt_diann.py` and `tmt_spectronaut.py`
> are in the same directory when executing!

### Chimerys DIA

- Export Chimerys PSMs from Proteome Discoverer in tab-separated .txt format.
- Set you desired parameters in `config.toml`
- Download ThermoRawFileParser from https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v1.4.5
- Convert your RAW file with `ThermoRawFileParser.exe -i RAW_FILE_NAME`
- Option A (recommended): Run via [uv](https://docs.astral.sh/uv/)
  - [Install uv](https://docs.astral.sh/uv/getting-started/installation/) if it's not already installed on your system, e.g.: `pip install uv`
  - Run the script with `uv run tmt_chimerys.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMS.txt -c config.toml`
- Option B: Run via native python
  - Install python 3.12 or greater from https://www.python.org/downloads/
  - Install requirements with `pip install -r requirements.txt`
  - Run the script with `python tmt_chimerys.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMS.txt -c config.toml`

### Chimerys DDA

- Export Chimerys PSMs from Proteome Discoverer in tab-separated .txt format.
- Set you desired parameters in `config.toml`
- Download ThermoRawFileParser from https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v1.4.5
- Convert your RAW file with `ThermoRawFileParser.exe -i RAW_FILE_NAME`
- Option A (recommended): Run via [uv](https://docs.astral.sh/uv/)
  - [Install uv](https://docs.astral.sh/uv/getting-started/installation/) if it's not already installed on your system, e.g.: `pip install uv`
  - Run the script with `uv run tmt_chimerys_dda.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMS.txt -c config.toml`
- Option B: Run via native python
  - Install python 3.12 or greater from https://www.python.org/downloads/
  - Install requirements with `pip install -r requirements.txt`
  - Run the script with `python tmt_chimerys_dda.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMS.txt -c config.toml`
