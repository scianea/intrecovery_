# PICRUSt2 functional profiling and RDA of anaerobic enrichment reactors

R code for two analyses of six anaerobic enrichment reactors (R1–R6). The microbial communities were profiled with Oxford Nanopore 16S rRNA sequencing and EMU, and their functions were predicted with PICRUSt2.

| Analysis | Folder | What it does |
|----------|--------|--------------|
| **PICRUSt2 analysis** | [`picrust2_analysis/`](picrust2_analysis/) | Differential abundance of MetaCyc pathways (DESeq2), run on all MetaCyc pathways and on an organonitrogen subset, plus a combined Z-score figure for the manuscript |
| **RDA** | [`rda_analysis/`](rda_analysis/) | Redundancy analysis linking treatment to process chemistry, community composition (EMU PCoA) and predicted function (PICRUSt2 PCA) |

Each folder has its own README with the full method and the outputs.

## Treatments

| Code | Treatment | Notes |
|------|-----------|-------|
| R1 | Control | |
| R2 | BES | Reference group for R3–R6 in the DESeq2 comparisons |
| R3 | Cysteine | |
| R4 | Ascorbate | |
| R5 | Titanium | |
| R6 | Mix | |

Treatment colours are fixed in `TREATMENT_COLORS` so R1–R6 look the same in every figure.

## Repository layout

```
.
├── picrust2_analysis/
│   ├── 01_picrust2_dual_analysis.R    # DESeq2 DAA, PCA, functional alpha diversity, heat maps
│   ├── 02_zscore_combined_figure.R    # combined Z-score figure (A/B) for the manuscript
│   └── README.md
├── rda_analysis/
│   ├── rda_analysis.R                 # RDA on both PICRUSt2 datasets + comparison table
│   └── README.md
├── install_packages.R                 # installs every required package
└── picrust2-rda-analysis.Rproj        # opens the repo with the right working directory
```

## Requirements

- R ≥ 4.2
- CRAN: `tidyverse`, `ggpicrust2`, `vegan`, `rstatix`, `ggpubr`, `tidytext`, `ragg`, `patchwork`, `readxl`, `ggrepel`
- Bioconductor: `DESeq2` (used internally by ggpicrust2)
- The **Times New Roman** font installed on the system (it is the default on Windows and macOS; on Linux install `ttf-mscorefonts-installer` or change `PUB_FONT` at the top of each script)

```r
source("install_packages.R")
```

## Quick start

All paths are relative to the repository root. Open `picrust2-rda-analysis.Rproj` in RStudio (or `setwd()` to the repo folder), put the input files in `data/` as described under [Input data](#input-data), then run:

```r
# 1. PICRUSt2 analysis
source("picrust2_analysis/01_picrust2_dual_analysis.R")
source("picrust2_analysis/02_zscore_combined_figure.R")

# 2. RDA (uses the PICRUSt2 PCA scores, so run it after the PICRUSt2 analysis)
source("rda_analysis/rda_analysis.R")
```

## Input data

This repository contains code only. To run it, create a `data/` folder in the repository root with this structure. Outputs are written to `results/`, which the scripts create automatically. Both folders are listed in `.gitignore`, so they are never uploaded.

```
data/
├── picrust2_output/                      # PICRUSt2 output folder
│   └── pathways_out/
│       └── path_abun_unstrat.tsv.gz      # used by the PICRUSt2 analysis
├── chemistry/
│   └── All_Results_RDA_M1.xlsx           # used by the RDA
└── emu_results/                          # used by the RDA
    ├── final_analysis/
    │   ├── selected_best_runs.csv        # columns: treatment, day, run, barcode
    │   └── PCoA_coordinates.csv          # columns: treatment, day, PCoA1, PCoA2
    ├── Batch1_Final/                     # EMU files for run "Batch1Final"
    ├── TiCtMix/                          # EMU files for run "TiCtMix"
    └── *_merged.fastq_rel-abundance.tsv  # EMU files for runs "CysAsc" and "CtrlBES"
```

### File details

**`path_abun_unstrat.tsv.gz`**: first column = MetaCyc pathway ID, one column per sample. Sample names must start with the treatment (`Control`, `BES`, `Cysteine`, `Ascorbic`, `Titanium`, `Mixed`) and contain the sampling day as the first number. Columns `SeedS1` and `SeedS2` are dropped automatically.

**`All_Results_RDA_M1.xlsx`**: four sheets, each with a `Treatment` column:

| Sheet | Day column | Variables used |
|-------|-----------|----------------|
| `NH4_NH3` | `Retention Time (day)` | `NH4_mean (mmol/L)`, `NH3_mean (mmol/L)` |
| `PH_ORP` | `Retention Time (day)` | `pH`, `ORP (mV)` (decimal commas are accepted) |
| `AA_Amines` | `Day` | `Total Amine Production (mM)`, `Total Amino Acid Production (mM)`, `Total Organonitrogen (mM)` |
| `SCFAs_Total C` | `Day` | `Total SCFA Concentration (mM)` |

**EMU files**: `<run>_<barcode>_merged.fastq_rel-abundance.tsv` with the columns `species` and `estimated counts`. The folder for each run is set in `run_to_folder` in `rda_analysis/rda_analysis.R`.

Treatment names in every file are standardised by `standardize_treatment_col()`. Both long names and the R1–R6 codes are accepted.

## Upstream steps (not in this repository)

- ONT 16S rRNA basecalling and **EMU** taxonomic profiling
- `emu_final_analysis` script, which selects the best run per sample (`selected_best_runs.csv`) and computes the community PCoA (`PCoA_coordinates.csv`)
- `picrust2_pipeline.py` (PICRUSt2), run on the EMU-derived abundance table of the 59 selected samples

## Citation

If you use this code, please cite the associated manuscript (reference to be added on publication).

## Contact

Ege Ulusman – University of Eastern Finland – eulusman@uef.fi
