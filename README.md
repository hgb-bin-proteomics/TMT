# TMT

TMTpro-18plex quantification for \[single cell\] DIA and DDA searches with
[Chimerys](https://www.msaid.de/chimerys),
[Spectronaut](https://biognosys.com/software/spectronaut/), and
[DIA-NN](https://github.com/vdemichev/DiaNN).

## Requirements

- Please install [OpenMS](https://openms.readthedocs.io/en/latest/about/installation.html).
  - We recommend and tested using OpenMS version [3.4.0](https://abibuilder.cs.uni-tuebingen.de/archive/openms/OpenMSInstaller/release/3.4.0/)!
- If you want to run the python scripts, you need to install [python 3.12 or higher](https://www.python.org/downloads/)
  or [uv](https://docs.astral.sh/uv/).

## Usage

- On Microsoft Windows the applications can be run as standalone executables or as python scripts.
- Other operating systems are limited to the python scripts, please refer to [CLI.md](CLI.md).

### Graphical User Interface

![GUI screenshot](docs/gui.png)

We provide compiled binaries for Microsoft Windows that offer a graphical user interface. Please download the executables from
[releases](https://github.com/hgb-bin-proteomics/TMT/releases).

> [!IMPORTANT]
>
> Please make sure that the executable and the `tmt18plex_default.ini` file are in the same directory!
> You might also have to unblock the `tmt18plex_default.ini` file either via its _Properties_ (right-click) or
> using [PowerShell](https://learn.microsoft.com/en-us/powershell/module/microsoft.powershell.utility/unblock-file).

### Commandline Interface

Please refer to [CLI.md](CLI.md).

### Running the Scripts for Multiple Files

If you want to run the scripts for multiple input files sequentially, please
refer to [MULTI.md](MULTI.md).

### Configuration

Please set the following parameters according to your needs in the `config.toml` file:

```toml
[METHOD]
# window size
window_size = 0.5
# window start (m/z)
window_start = 400.0
# window end (m/z)
window_end = 800.0
# window overlap
window_overlap = 0.0

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
# precursor intensity fraction in the window to use as reference, only used for displaying some preliminary statistics
# for filtering use PROTEIN.min_purity
total_intensity_threshold = 0.7
# minimum relative intensity threshold compared to most intense peak in window to not be considered noise
noise_threshold = 0.1

[QUANTIFICATION]
# subtract the reporter noise from the reporter signal?
# if true, filtering by S/N should be turned off or thresholds set to 0.0
subtract_noise = true
# quantification method to use
# 1 = native
# 2 = OpenMS
# 3 = Resolution GUI
quantification_method = 2

[PROTEIN]
# Qvalue that should be used for filtering, only used for DIA-NN and Spectronaut
q_value = 0.01
# Normalized Chimerys Coefficient Threshold, anything below will be ignored, only applies to Chimerys
min_chimerys_coefficient = 1.0
# minimum average reporter S/N for a PSM to be considered for aggregation, only applies to Chimerys
min_avg_reporter_sn = 10.0
# minimum reporter resolution to be considered for aggregation
min_reporter_res = 45000.0
# minimum purity for a PSM to be considered for aggregation
min_purity = 0.7
# whether or not ambiguous protein groups should be filtered out, only used for DIA-NN and Spectronaut
keep_ambiguous_protein_groups = false

[CONDITIONS]
# please define your conditions here
# conditions should be given as condition name (without spaces) following an equal sign and then a list of TMT reporters
# see examples below
all = ["TMTpro-126",  "TMTpro-127N", "TMTpro-127C", "TMTpro-128N", "TMTpro-128C",
       "TMTpro-129N", "TMTpro-129C", "TMTpro-130N", "TMTpro-130C", "TMTpro-131N",
       "TMTpro-131C", "TMTpro-132N", "TMTpro-132C", "TMTpro-133N", "TMTpro-133C",
       "TMTpro-134N", "TMTpro-134C", "TMTpro-135N"]
cond1 = ["TMTpro-126",  "TMTpro-127N", "TMTpro-127C", "TMTpro-128N", "TMTpro-128C",
         "TMTpro-129N", "TMTpro-129C", "TMTpro-130N", "TMTpro-130C"]
cond2 = ["TMTpro-131N", "TMTpro-131C", "TMTpro-132N", "TMTpro-132C", "TMTpro-133N",
         "TMTpro-133C", "TMTpro-134N", "TMTpro-134C", "TMTpro-135N"]
# please define the min S/N thresholds per condition that should be used for protein aggregation here
# this should be sn_thresholds = map of thresholds for each condition
# see example below
sn_thresholds = { all = 0.0, cond1 = 10.0, cond2 = 10.0 }
```

> [!IMPORTANT]
>
> You might also want to adapt the isotope correction factors for your TMT lot, you can do that in the `tmt18plex_default.ini` file.
> Please refer to the documentation site of OpenMS [here](https://openms.de/documentation/html/TOPP_IsobaricAnalyzer.html).

## TMT Resolution GUI Tool

You might also want to use the output of the Resolution GUI tool developed by Dina L. Bai, Tian Zhang _et al._ [\[1\]](https://doi.org/10.1038/s41467-025-60022-x) as additional input for better quality control. Please refer to this repository for instructions: [https://github.com/hgb-bin-proteomics/TMT_Resolution_GUI](https://github.com/hgb-bin-proteomics/TMT_Resolution_GUI).

## Contact

In case of questions please contact:
- [micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at)
