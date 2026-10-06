# PICRUSt2 analysis

Differential abundance analysis of PICRUSt2-predicted MetaCyc pathways, run on two datasets:

1. **All MetaCyc** pathways
2. **Organonitrogen** pathways: a subset selected by keyword match on the MetaCyc descriptions (amino acids, amines and polyamines, urea, heme/tetrapyrroles, nitrogen-cycle terms)

The analysis has two scripts that run in order. Both are run from the repository root.

| Step | Script | Input | Output folder |
|------|--------|-------|---------------|
| 1 | `01_picrust2_dual_analysis.R` | `data/picrust2_output/pathways_out/path_abun_unstrat.tsv.gz` | `results/picrust2_analysis/{All_MetaCyc,Organonitrogen}/` |
| 2 | `02_zscore_combined_figure.R` | step 1 results (`all_daa_results.csv`) + the same abundance file | `results/picrust2_analysis/Combined_Figures/` |

```r
source("picrust2_analysis/01_picrust2_dual_analysis.R")
source("picrust2_analysis/02_zscore_combined_figure.R")
```

## Step 1: differential abundance and exploratory plots

For each dataset:

- **DESeq2** (through `ggpicrust2::pathway_daa`) pairwise comparisons R1→R2, R2→R3, R2→R4, R2→R5, R2→R6; significance at adjusted p < 0.05. Pathways with zero counts in a comparison are removed before testing.
- **PCA** on centred and scaled pathway abundances. The PC1/PC2 scores are saved and used later by the RDA.
- **Functional alpha diversity** (richness, Shannon, Pielou): Kruskal–Wallis test, then BH-adjusted pairwise Wilcoxon tests.
- **Figures:** log2 fold-change heat map (top 30 pathways), top-15 bar plot per comparison, per-treatment Z-score heat map (top 15 per treatment).

Sample names are parsed into treatment and day; the script stops with an error if a sample cannot be matched. Seed samples (`SeedS1`, `SeedS2`) are excluded.

| Output | Content |
|--------|---------|
| `all_daa_results.csv` | All DESeq2 results with MetaCyc descriptions |
| `significant_pathways.csv` | Hits with adjusted p < 0.05 |
| `PCA_PC1_PC2.csv` | PCA scores per sample (input to the RDA) |
| `Alpha_Diversity.csv`, `Alpha_Kruskal_Wallis.csv`, `Alpha_Pairwise_Wilcoxon.csv` | Alpha diversity values and tests |
| `PCA.png`, `Functional_Alpha_Diversity.png`, `Shannon_Diversity.png`, `fold_change_heatmap.png`, `Top15_Pathways_Per_Comparison.png`, `zscore_heatmap.png` | Figures (300 dpi) |
| `Organonitrogen/organonitrogen_pathways.txt` | The pathways in the organonitrogen subset |

## Step 2: combined Z-score figure (manuscript figure)

- Reuses `all_daa_results.csv` from step 1, so DESeq2 does not run again.
- Takes the top 10 significant pathways per treatment and shows within-treatment Z-scores across sampling days.
- Panel **A** = All MetaCyc, panel **B** = Organonitrogen, with one shared, symmetric colour scale.
- Labels are kept on one line, HTML entities become Greek letters, and some long names are abbreviated. Define these abbreviations in the figure legend: SP, superpathway of; GlcNAc, N-acetylglucosamine; ManNAc, N-acetylmannosamine; Neu5Ac, N-acetylneuraminate.
- Saved at 7.5 × 9.5 in (full page width) as PNG and LZW-compressed TIFF at 600 dpi plus a vector PDF.

| Output | Content |
|--------|---------|
| `Figure_ZScore_Combined.{png,tiff,pdf}` | The figure |
| `Figure_ZScore_Combined_data.csv` | Plotted data, including the full pathway names behind each label |

## Settings you may want to change

At the top of each script:

- `ALPHA_CUTOFF` is the significance threshold (default 0.05).
- `TOP_N_HEATMAP`, `TOP_N_BARS` and `TOP_N_ZSCORE` set the number of pathways shown.
- `comparisons` lists the DESeq2 contrasts.
- `FIG_WIDTH`/`FIG_HEIGHT` (step 2) set the printed size; check your journal's limits.
- `LABEL_MODE` (step 2): `"name"` shows pathway names, `"id"` shows MetaCyc IDs.
