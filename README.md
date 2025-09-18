# TMT

TMTpro-18plex quantification for \[single cell\] DIA and DDA searches with
[Chimerys](https://www.msaid.de/chimerys),
[Spectronaut](https://biognosys.com/software/spectronaut/), and
[DIA-NN](https://github.com/vdemichev/DiaNN).

## Usage

- On Microsoft Windows the applications can be run as standalone executables or as python scripts.
- Other operating systems are limited to the python scripts, please refer to [CLI.md](CLI.md).

### Graphical User Interface

![GUI screenshot](docs/gui.png)

We provide compiled binaries for Microsoft Windows that offer a graphical user interface. Please download the executables from
either [releases](https://github.com/hgb-bin-proteomics/TMT/releases) or in zipped form from the
[binaries](https://github.com/hgb-bin-proteomics/TMT/tree/master/binaries) folder.

> [!IMPORTANT]
>
> Please make sure that the executable and the `tmt18plex_default.ini` file are in the same directory.

### Commandline Interface

Please refer to [CLI.md](CLI.md).

### Configuration

Please set the following parameters according to your needs in the `config.toml` file.

```toml
[METHOD]
# window size
window_size = 0.5
# window start (m/z)
window_start = 400.0
# window end (m/z)
window_end = 800.0

[MATCHING]
# m/z tolerance for matching peaks in Dalton
mz_tolerance = 0.02
# retention time tolerance in seconds for matching identifications to MS2 spectra
rt_tolerance = 3.0
# retention time window in seconds that a MS1 and corresponding MS2 spectrum must be in
ms1_rt_window = 10.0

[ISOTOPES]
# whether precursor isotopes should be considered for purity calculation
consider_precursor_isotopes = true
# isotope match tolerance in Dalton
isotope_tolerance = 0.01
# maximum considered precursor charge
max_charge = 6

[FILTERING]
# precursor intensity fraction in the window to use as reference
total_intensity_threshold = 0.7
# minimum relative intensity threshold compared to most intense peak in window to not be considered noise
noise_threshold = 0.1

[PROTEIN]
# Normalized Chimerys Coefficient Threshold, anything below will be ignored
min_chimerys_coefficient = 1.0
# minimum average reporter S/N for a PSM to be considered for aggregation
min_avg_reporter_sn = 10.0
# minimum reporter resolution to be considered for aggregation
min_reporter_res = 45000.0
# minimum purity for a PSM to be considered for aggregation
min_purity = 0.7
```

## Contact

In case of questions please contact:
- [micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at)
