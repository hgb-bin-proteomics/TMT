# TMT

## Usage

- Download ThermoRawFileParser from https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v1.4.5
- Convert your RAW file with `ThermoRawFileParser.exe -i RAW_FILE_NAME`
- Install python 3.12 or greater from https://www.python.org/downloads/
- Install requirements with `pip install -r requirements.txt`
- Set you desired parameters in `config.toml`
- Export Chimerys PSMs from Proteome Discoverer in Microsoft Excel format.
- Run the script with `python tmt_chimerys.py -s SPECTRA.mzML -i PROTEOME_DISCOVERER_PSMS.xlsx -c config.toml`
