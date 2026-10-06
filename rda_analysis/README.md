# RDA

Redundancy analysis (RDA, `vegan::rda`) that tests how much of the variation in process chemistry, community composition and predicted function is explained by treatment.

**Run the [PICRUSt2 analysis](../picrust2_analysis/) first.** This script reads its `PCA_PC1_PC2.csv` files.

```r
source("rda_analysis/rda_analysis.R")
```

## Model

- **Explanatory variable:** Treatment (R1–R6)
- **Response variables** (all scaled to unit variance):

| Group | Variables | Source |
|-------|-----------|--------|
| Chemistry | pH, ORP, total SCFA, TAN (NH₄⁺ + NH₃), FAN (NH₃), total amino acids, total amines | `data/chemistry/All_Results_RDA_M1.xlsx` |
| Community | EMU PCoA1, PCoA2 | `data/emu_results/final_analysis/PCoA_coordinates.csv` |
| Function | PICRUSt2 PC1, PC2 | `results/picrust2_analysis/<dataset>/PCA_PC1_PC2.csv` |

The model is fitted twice, once with the All MetaCyc PCs and once with the Organonitrogen PCs. Everything else is identical, so the two models can be compared directly.

- Only complete cases are used: 4 of the 59 samples lack NH₄⁺/NH₃ data and are dropped.
- Significance comes from permutation tests (`anova.cca`, 999 permutations, `set.seed(123)`) for the whole model, by term and by axis.
- Treatment names from every input are standardised in one place (`standardize_treatment_col()`), and the script checks that chemistry and community samples line up.

## Outputs

Written to `results/rda_analysis/`, with one subfolder per dataset (`All_MetaCyc/`, `Organonitrogen/`):

| Output | Content |
|--------|---------|
| `RDA_summary_and_tests.txt` | Model summary and permutation tests (overall, by term, by axis) |
| `RDA_input_data.csv` | The exact data the model was fitted on |
| `rda1v2.png`, `rda1v3.png` | Triplots: samples with 95% ellipses, treatment centroids, response arrows; axis labels give % variance and the p-value per axis |
| `rda_triplot_combined.png` | RDA1 vs RDA2 (A) above RDA1 vs RDA3 (B), with a shared legend |

In the `results/rda_analysis/` root:

| Output | Content |
|--------|---------|
| `RDA_comparison_All_vs_Organonitrogen.csv` | Constrained variance (%), adjusted R², pseudo-F and p for both models |

## Note

The script also builds the full EMU species count matrix. It is used only to check that the samples line up and is kept for possible future analyses. The model itself uses PCoA1/PCoA2 as the community summary.
