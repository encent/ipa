# ipa

**Interaction Pattern Aggregation analysis**

The goal of this analysis is to assess the average pattern shape (e.g. loops or micro-compartments) in the 3C data, and, if necessary, to compare it with the other signals: ATAC-Seq, histone modifications, etc. 

***ipa* is useful when:**

1. You want to assess the shape of a 3C pattern in addition to the classic pileup analysis (using ***cooltools*** or ***coolpup.py***). For example, to show the difference in the shape of loops and micro-compartments (Fig. B).
2. You want to compare the shape of a 3C pattern (e.g. loops or micro-compartments) with some other signals, e.g. ATAC-Seq, histone modifications (Fig. C). 

**The analysis is done in two steps:**

1. The sum of contacts is calculated for every bin in the ICE-normalized or in the observed-over-expected 3C matrix. Calculation of the sum of contacts is restricted to the minimum and maximum size of a pattern which shape should be assessed, e.g. minimum/maximum loop size (Fig. A). The output of this step is a bigWig track for the calculated restricted sum of contacts.
2. Using [***pybbi***](https://github.com/nvictus/pybbi), the stackup plots (Fig. B—C) are created based on the input bigWig files and the annotation file with the genomic coordinates of the loci that are located in a close proximity to the examined pattern: loop anchors, domain borders, TSS/TES, etc. If the annotation file has the strand column, the corresponding intervals from the bigWig file will be properly flipped.

**Installation and usage**

To install ***ipa*** on your machine, follow the instructions from the sections [System Requirements](#system-requirements) and [Getting Started](#getting-started).

If you want to reproduce the ***ipa*** analysis from the *Kim et al.* paper, go to the section [Example: reproducing ***ipa*** plots from the *Kim et al.* paper](#example-reproducing-ipa-plots-from-the-kim-et-al-paper).

If you want to run ***ipa*** on your own datasets, please go to the section [Usage](#usage).

![overview](images/overview.png)

## System Requirements

You must have Python 3 with `pip` and `conda` installed on your system:
```bash
python3 --version
python3 -m pip --version
conda --version
```
If Python 3 is missing, follow the installation instructions from [the Python official website](https://www.python.org/). If you have Python 3.4 or later, `pip` is included by default.

If `conda` is missing, follow [these instructions](https://www.anaconda.com/docs/getting-started/miniconda/install) to install Miniconda.

## Getting Started

To get started with ***ipa***, follow these steps:

1. **Clone the repository**:
  ```bash
  git clone https://github.com/encent/ipa.git
  cd ipa
  ```

2. **Create conda environment**:
  ```bash
  conda create -n ipa python=3.12
  ```

3. **Activate conda environment**:
  ```bash
  conda activate ipa
  ```

4. **Install the dependencies**:
  ```bash
  pip install -r requirements.txt
  ```

5. **Install *ipa***:
  ```bash
  pip install .
  ```

## Usage

### Command-Line Interface

***ipa*** has three commands:

#### `ipa track`

This command calculates the sum of contacts, resticted to the minimum/maximum pattern size, for every bin in the ICE-normalized or in the observed-over-expected 3C matrix, and produces a corresponding bigWig file.

**Usage:**

```bash
ipa track [OPTIONS]
```

**Options:**

* `--cool-path`, `--cool_path`, `-c` **(required)**:
                        Path to the .cool file. Keep in mind that chromosome names in the .cool file and in all the files that will be used in the `ipa plot` (e.g. TSS-TES sites, ATAC-Seq signal bigWig file) should match each other.
* `--output-dir`, `--output_dir`, `-o` **(required)**:
                        Path to create the output directory which will store the output bigWig file.
* `--expected`, `-e`: 
                        If `True`, generates an ***ipa*** track based on the observed over expected matrix (default: `False`).
* `--clr-weight-name`, `--clr_weight_name`, `-b`:
                        The name of the column in the cool file that contains the balancing weights (default: `'weight'`).
* `--min-dist`, `--min_dist`:
                        Minimum distance (in bp) between two loci to consider for the ***ipa*** calculation, e.g. minimum loop size in bp. If `None`, restriction on minimum distance is not applied (default: `40_000`).
* `--max-dist`, `--max_dist`:
                        Maximum distance (in bp) between two loci to consider for the ***ipa*** calculation, e.g. maximum loop size in bp. If `None`, restriction on maximum distance is not applied (default: `100_000`).
* `--nproc`, `-np`:
                        Number of processes to use for the calculation of expected. Used when `--expected` is `True` (default: `4`).

**Example:**

```bash
ipa track \
          --cool-path /path/to/cool/file.mcool::resolutions/5000 \
          --output-dir /path/to/output/dir \
          --clr-weight-name weight \
          --min-dist 40000 \
          --max-dist 1000000 \
          --nproc 4
```

#### `ipa plot`

This command creates a stackup plot from one or two bigWig files (including one that was generated from the 3C matrix by the `ipa track` command) and the set of regions of interest: genes, domains, loop coordinates etc.

**Usage:**

```bash
ipa plot [OPTIONS]
```

**Options:**
* `--bw-path`, `--bw_path`, `-bw` **(required)**:
                        Path to the bigWig file. Keep in mind that chromosome names in the .bw files and in all the files that will be used in the ipa_plot() function (e.g. TSS-TES sites, ATAC-Seq signal .bw file) should match each other.
* `--roi-path`, `--roi_path`, `-roi` **(required)**:
                        Path to the annotation file with the regions of interest (e.g. TSS-TES sites). The file should be in a BED format.
* `--output-dir`, `--output_dir`, `-o` **(required)**:
                        Path to the output directory which will store the output plot file.
* `--extra-bw-path`, `--extra_bw_path`, `-extra`:
                        Path to the second bigWig file (default: `None`).
* `--roi-start-name`, `--roi_start_name`:
                        Alias for the start of the region of interest, e.g. TSS or loop start (default: `None`).
* `--roi-end-name`, `--roi_end_name`:
                        Alias for the end of the region of interest, e.g. TES or loop end (default: `None`).
* `--flank`:            Size of the flanking regions in bp (default: `100_000`).
* `--nbins`:            Number of bins for the stackup plot (default: `50`).
* `--min-roi-size`, `--min_roi_size`:
                        Minimum size of the region of interest (ROI) in bp to filter out small regions in the roi file (default: `None`).
* `--max-roi-size`, `--max_roi_size`:
                        Maximum size of the region of interest (ROI) in bp to filter out large regions in the roi file (default: `None`).

**Example:**

```bash
ipa plot \
          --bw-path /path/to/sum/of/contacts/file.bw \
          --roi-path /path/to/roi/file.bed \
          --output-dir /path/to/output/dir \
          --extra-bw-path /path/to/extra.bw \
          --roi-start-name TSS \
          --roi-end-name TES \
          --flank 100000 \
          --nbins 50
```

#### `ipa`

This is the main command that runs the entire analysis. It first runs the `ipa track` command and then `ipa plot` command. The main advantage of running this command, compared to running two other commands consequtively, is that it can plot as many stackups as there are additional bigWig files in the corresponding folder. This is useful when you either have many species/experiments to run in a row for different sets of parameters, or when you many bigWig files which you want to plot in addition to the 3C-based stackup. 

**Usage:**

```bash
ipa [OPTIONS]
```

**Options:**
* `--cool-path`, `--cool_path`, `-c` **(required)**:
                        Path to the .cool file. Keep in mind that chromosome names in the .cool file and in all the files that will be used in the `ipa plot` (e.g. TSS-TES sites, ATAC-Seq signal bigWig file) should match each other.
* `--roi-path`, `--roi_path`, `-roi` **(required)**:
                        Path to the annotation file with the regions of interest (e.g. TSS-TES sites). The file should be in a BED format.
* `--output-dir`, `--output_dir`, `-o` **(required)**:
                        Path to the output directory which will store the output bigWig file and the output plot files.
* `--bw-dir`, `--bw_dir`, `-bw`:
                        Path to the directory with the bigWig files (ATAC-Seq, ChIP-Seq etc.) that will be used in the `ipa plot` function (default: `None`).
* `--expected`, `-e`: 
                        If `True`, generates an ***ipa*** track based on the observed over expected matrix (default: `False`).
* `--clr-weight-name`, `--clr_weight_name`, `-b`:
                        The name of the column in the cool file that contains the balancing weights (default: `'weight'`).
* `--min-dist`, `--min_dist`:
                        Minimum distance (in bp) between two loci to consider for the ***ipa*** calculation, e.g. minimum loop size in bp. If `None`, restriction on minimum distance is not applied (default: `40_000`).
* `--max-dist`, `--max_dist`:
                        Maximum distance (in bp) between two loci to consider for the ***ipa*** calculation, e.g. maximum loop size in bp. If `None`, restriction on maximum distance is not applied (default: `100_000`).
* `--nproc`, `-np`:
                        Number of processes to use for the calculation of expected. Used when `--expected` is `True` (default: `4`).
* `--roi-start-name`, `--roi_start_name`:
                        Alias for the start of the region of interest, e.g. TSS or loop start (default: `None`).
* `--roi-end-name`, `--roi_end_name`:
                        Alias for the end of the region of interest, e.g. TES or loop end (default: `None`).
* `--flank`:            Size of the flanking regions in bp (default: `100_000`).
* `--nbins`:            Number of bins for the stackup plot (default: `50`).
* `--min-roi-size`, `--min_roi_size`:
                        Minimum size of the region of interest (ROI) in bp to filter out small regions in the roi file (default: `None`).
* `--max-roi-size`, `--max_roi_size`:
                        Maximum size of the region of interest (ROI) in bp to filter out large regions in the roi file (default: `None`).

**Example:**

```bash
ipa \
    --cool-path /path/to/cool/file.mcool::resolutions/5000 \
    --roi-path /path/to/roi/file.bed \
    --output-dir /path/to/output/dir \
    --bw-dir /path/to/dir/with/bigwig/files \
    --clr-weight-name weight \
    --min-dist 40000 \
    --max-dist 1000000 \
    --nproc 4 \
    --roi-start-name TSS \
    --roi-end-name TES \
    --flank 100000 \
    --nbins 50
```

### Python API

See the `ipa.py` docstrings for more details.

## Example: reproducing `ipa` plots from the *Kim et al.* paper

To reproduce the ***ipa*** plots from the paper, first install ***ipa*** ([System Requirements](#system-requirements) and [Getting Started](#getting-started)). After that, download the input data using [this link](https://ccnag-my.sharepoint.com/:f:/g/personal/nikolai_bykov_cnag_eu/ElwkqPgVuTdAoWlEUXUwICYBnYK74_WPbbipFTryXlJlQg?e=47sSwg). All datasets in this folder are ours except for those related to human (Krietenstein et al., 2020) and fruit fly (Batut et al., 2022). The full list of public datasets with accession numbers used in our study is provided in the Supplementary Table 2 in the *Kim et al.* paper.

Go to the `ipa/scripts` folder:

```bash
cd scripts
```

Run the `Kim_et_al_script.py` script:

```bash
python3 Kim_et_al_script.py \
        --input-dir /path/to/the/downloaded/data_Kim_et_al \
        --output-dir /path/to/the/output/directory
```

The `--input-dir` parameter should be a path to the folder with downloaded data. The plots will be placed in the path you provide in the parameter `--output-dir`. Their design will differ from those on the figure above, but all the data should be the same.

***NB:*** The script has been tested on a Linux machine with 128 GB of RAM. Running time was approximately 21 minutes (for the seven species and the parameters provided in the script). Please be informed that the current version of ***ipa*** could be memory-consuming for large genomes, such as human (the script is processing human Micro-C data). Thus, if you have less than 128 GB of RAM on your machine, we cannot guarantee that the script will run perfectly on human data.

### `ipa` parameters used in this script

Sum of contacts was calculated over the ICE-normalized matrix. For the ***ipa*** plots we used parameter `--nbins 50`, which indicates the number of bins to build a stackup plot using [***pybbi***](https://github.com/nvictus/pybbi).

For TSS-TES comparative plots (Fig. C) we used these species-specific ***ipa*** parameters:
| Species | Micro-C resolution (bp) | Size range of contacts (Kb) | Flank (Kb) |
|---------|-------------------------|-----------------------------|------------|
| Cowc    | 400                     | 4 — 100                     | 10         |
| Tadh    | 400                     | 4 — 100                     | 10         |
| Mlei    | 800                     | 5 — 150                     | 10         |
| Nvec    | 500                     | 10 — 360                    | 10         |
| Dmel    | 400                     | 5 — 250                     | 10         |
| Hsap    | 5000                    | 50 — 1060                   | 100        |

Size range of contacts (Kb) corresponds to the ***ipa*** parameters `--min-dist` and `--max-dist`. 


To obtain the domain pattern shape plot for *Sphaeroforma arctica* (Fig. B) we used these species-specific ***ipa*** parameters:
| Species | Micro-C resolution (bp) | Size range of contacts (Kb) | Flank (Kb) |
|---------|-------------------------|-----------------------------|------------|
| Sarc    | 2800                    | 50 — 4000                   | 30         |

## Documentation

Documentation is currently provided in the docstrings and in this README.

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue for any enhancements or bug fixes.

## Citing `ipa`

[TBD] Kim et al., 2025
