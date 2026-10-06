# ==============================================================================
# v2 - same analysis as picrust2_first_manuscript_edited_3009.R, re-run on the
# PICRUSt2 output rebuilt from the 59 selected EMU runs (picrust2_output_munia1_v2).
# Original script, input folder and results are kept unchanged for comparison.
# ==============================================================================
# PICRUSt2 -> ggpicrust2/DESeq2 differential abundance: DUAL ANALYSIS
#   1. All MetaCyc pathways
#   2. Organonitrogen MetaCyc pathways (subset)
# Run AFTER picrust2_pipeline.py has completed.
# ==============================================================================

library(tidyverse)   # ggplot2, dplyr, tidyr, readr, purrr, stringr, tibble
library(ggpicrust2)  # loads DESeq2 internally; no need to attach it
library(vegan)
library(rstatix)
library(ggpubr)
library(tidytext)
library(ragg)

# ------------------------------------------------------------------------------
# 1. Config
# ------------------------------------------------------------------------------

picrust2_dir <- "data/picrust2_output"              # PICRUSt2 output folder (contains pathways_out/)
results_root <- "results/picrust2_analysis"    # created if missing

# R1 = Control, R2 = BES (reference for R3-R6), R3 = Cysteine, R4 = Ascorbic,
# R5 = Titanium, R6 = Mixed
comparisons <- list(
  c("R1", "R2"),
  c("R2", "R3"),
  c("R2", "R4"),
  c("R2", "R5"),
  c("R2", "R6")
)
COMPARISON_LEVELS <- map_chr(comparisons, ~ paste0(.x[1], " → ", .x[2]))

ALPHA_CUTOFF <- 0.05
TOP_N_HEATMAP <- 30  # pathways in the fold-change heatmap
TOP_N_BARS    <- 15  # pathways per comparison in the bar plot
TOP_N_ZSCORE  <- 15  # pathways per treatment in the z-score heatmap

# ------------------------------------------------------------------------------
# 2. Publication theme (Times New Roman, one size scale)
# ------------------------------------------------------------------------------
# ragg resolves installed system fonts via systemfonts, so no registration is
# needed. ggsave() uses ragg automatically for .png when ragg is installed.

PUB_FONT      <- "Times New Roman"
PUB_BASE_SIZE <- 12

theme_pub <- function(base_size = PUB_BASE_SIZE, base_family = PUB_FONT) {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title    = element_text(family = base_family, face = "bold",
                                   size = base_size + 3, hjust = 0.5,
                                   margin = margin(b = 6)),
      plot.subtitle = element_text(family = base_family, size = base_size + 1,
                                   hjust = 0.5, margin = margin(b = 6)),
      plot.caption  = element_text(family = base_family, size = base_size - 3),
      axis.title    = element_text(family = base_family, size = base_size),
      axis.text     = element_text(family = base_family, size = base_size - 2),
      legend.title  = element_text(family = base_family, size = base_size),
      legend.text   = element_text(family = base_family, size = base_size - 1),
      strip.text    = element_text(family = base_family, size = base_size,
                                   face = "bold")
    )
}
theme_set(theme_pub())
update_geom_defaults("text", list(family = PUB_FONT))  # geom_text/annotations

# Locked treatment colours: identical to TREATMENT_COLORS in
# emu_final_analysis.R, so R1-R6 look the same in every manuscript figure
# (PCoA, PCA, alpha diversity, RDA) regardless of which levels are present.
TREATMENT_COLORS <- c(
  "R1" = "#F8766D",  # Control
  "R2" = "#B79F00",  # BES
  "R3" = "#00BA38",  # Cysteine
  "R4" = "#00BFC4",  # Ascorbate
  "R5" = "#619CFF",  # Titanium
  "R6" = "#F564E3"   # Mix
)

# Shared diverging scale for every fold-change / z-score tile or bar
scale_fill_diverging <- function(...) {
  scale_fill_gradient2(low = "red", mid = "white", high = "darkblue",
                       midpoint = 0, ...)
}

# ------------------------------------------------------------------------------
# 3. Load data, metadata, organonitrogen subset
# ------------------------------------------------------------------------------

parse_sample_name <- function(sample_name) {
  tibble(sample = sample_name) %>%
    mutate(
      treatment = case_when(
        str_starts(sample, "Control")  ~ "R1",
        str_starts(sample, "BES")      ~ "R2",
        str_starts(sample, "Cysteine") ~ "R3",
        str_starts(sample, "Ascorbic") ~ "R4",
        str_starts(sample, "Titanium") ~ "R5",
        str_starts(sample, "Mixed")    ~ "R6",
        TRUE ~ NA_character_
      ),
      day = readr::parse_number(sample)
    )
}

pathway_abundance <- read_tsv(
  file.path(picrust2_dir, "pathways_out/path_abun_unstrat.tsv.gz"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("SeedS1", "SeedS2")))

id_col <- colnames(pathway_abundance)[1]

metadata <- parse_sample_name(colnames(pathway_abundance)[-1]) %>%
  mutate(treatment = factor(treatment, levels = paste0("R", 1:6)))

# Fail early if a column name did not map to a treatment/day
bad_samples <- metadata %>% filter(is.na(treatment) | is.na(day)) %>% pull(sample)
if (length(bad_samples) > 0) {
  stop("Samples not matched by parse_sample_name(): ",
       paste(bad_samples, collapse = ", "))
}
print(table(metadata$treatment))

# Organonitrogen subset. Terms already covered by a broader substring were
# removed (e.g. "amine" covers ethylamine/propylamine/polyamine...,
# "nitrogen" covers nitrogen fixation/assimilation/metabolism), so the set of
# matched pathways is identical to the long list.
organonitrogen_terms <- c(
  # amino acids
  "amino acid", "arginine", "ornithine", "glutamine", "glutamate", "lysine",
  "methionine", "cysteine", "histidine", "tryptophan", "alanine", "tyrosine",
  "leucine", "valine", "threonine", "serine", "glycine", "aspartate",
  "asparagine", "proline",
  # amines & polyamines
  "amine", "cadaverine", "agmatine", "putrescine", "spermidine", "spermine",
  # other N-containing metabolites
  "urea", "heme", "tetrapyrrole", "porphyrin",
  # nitrogen cycling
  "nitrogen", "nitrification", "denitrification", "anammox", "nitrate",
  "nitrite", "ammonia", "ammonium", "nitric oxide", "nitrous oxide"
)

pathway_info <- pathway_annotation(
  pathway        = "MetaCyc",
  daa_results_df = data.frame(feature = pathway_abundance[[id_col]]),
  ko_to_kegg     = FALSE
)

organonitrogen_pathways <- pathway_info %>%
  filter(grepl(paste(organonitrogen_terms, collapse = "|"),
               description, ignore.case = TRUE))

cat("Organonitrogen pathways matched:", nrow(organonitrogen_pathways), "of",
    nrow(pathway_abundance), "\n")

organonitrogen_abundance <- pathway_abundance %>%
  filter(.data[[id_col]] %in% organonitrogen_pathways$feature)

# ------------------------------------------------------------------------------
# 4. Helpers (one job each)
# ------------------------------------------------------------------------------

save_plot <- function(plot, dir, name, width, height) {
  ggsave(file.path(dir, name), plot, width = width, height = height, dpi = 300)
}

# ---- Differential abundance ---------------------------------------------------

run_pairwise_daa <- function(abundance_matrix, metadata, ref, comp) {
  sub_meta <- metadata %>%
    filter(treatment %in% c(ref, comp)) %>%
    mutate(treatment = relevel(factor(as.character(treatment)), ref = ref))

  sub_abun <- abundance_matrix[, sub_meta$sample, drop = FALSE]
  sub_abun <- sub_abun[rowSums(sub_abun) > 0, , drop = FALSE]  # DESeq2 needs signal

  pathway_daa(
    abundance  = sub_abun,
    metadata   = sub_meta,
    group      = "treatment",
    daa_method = "DESeq2"
  )
}

run_all_daa <- function(abundance_matrix, metadata, comparisons) {
  results <- map(comparisons, function(cmp) {
    tryCatch(
      run_pairwise_daa(abundance_matrix, metadata, cmp[1], cmp[2]),
      error = function(e) {
        warning("DAA failed for ", cmp[1], " vs ", cmp[2], ": ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
  })
  results <- bind_rows(results)
  if (nrow(results) == 0) stop("No DAA results were produced.")

  pathway_annotation(pathway = "MetaCyc", daa_results_df = results,
                     ko_to_kegg = FALSE)
}

# Mean log2FC per pathway x comparison, significant hits only
prep_fold_changes <- function(daa_annotated) {
  daa_annotated %>%
    filter(p_adjust < ALPHA_CUTOFF, !is.na(description)) %>%
    group_by(feature, description, group1, group2) %>%
    summarise(log2_fc = mean(log2_fold_change), .groups = "drop") %>%
    mutate(comparison = factor(paste0(group1, " → ", group2),
                               levels = COMPARISON_LEVELS))
}

# ---- PCA ------------------------------------------------------------------------

plot_pca <- function(abundance_matrix, metadata, results_dir, analysis_label) {
  keep    <- apply(abundance_matrix, 1, var, na.rm = TRUE) > 0
  pca_res <- prcomp(t(abundance_matrix[keep, , drop = FALSE]),
                    center = TRUE, scale. = TRUE)
  var_exp <- 100 * pca_res$sdev^2 / sum(pca_res$sdev^2)

  scores <- as.data.frame(pca_res$x[, 1:2]) %>%
    rownames_to_column("sample") %>%
    left_join(metadata, by = "sample")
  write_csv(scores, file.path(results_dir, "PCA_PC1_PC2.csv"))

  # Built exactly like the PCoA in emu_final_analysis.R (points size 4,
  # filled 95% ellipses at alpha 0.2 drawn after the points)
  p <- ggplot(scores, aes(PC1, PC2, color = treatment)) +
    geom_point(size = 4) +
    stat_ellipse(aes(fill = treatment), geom = "polygon",
                 alpha = 0.2, color = NA) +
    scale_color_manual(values = TREATMENT_COLORS) +
    scale_fill_manual(values = TREATMENT_COLORS) +
    labs(
      title    = paste(analysis_label, "- PCA of Functional Profiles"),
      subtitle = "PICRUSt2 MetaCyc pathway abundances grouped by treatment",
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2]),
      color = "Treatment", fill = "Treatment"
    ) +
    theme_pub()

  save_plot(p, results_dir, "PCA.png", 8, 6)
}

# ---- Alpha diversity + statistics ----------------------------------------------

run_alpha_diversity <- function(abundance_matrix, metadata, results_dir,
                                analysis_label) {
  m <- as.matrix(abundance_matrix)

  alpha_df <- tibble(
    sample   = colnames(m),
    Richness = colSums(m > 0),
    Shannon  = as.numeric(vegan::diversity(t(m), index = "shannon"))
  ) %>%
    mutate(Pielou = if_else(Richness > 1, Shannon / log(Richness), NA_real_)) %>%
    left_join(metadata, by = "sample")
  write_csv(alpha_df, file.path(results_dir, "Alpha_Diversity.csv"))

  metrics <- c("Richness", "Shannon", "Pielou")
  alpha_long <- alpha_df %>%
    pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = metrics))

  # Kruskal-Wallis per metric, then BH-adjusted pairwise Wilcoxon
  kw_results <- alpha_long %>%
    group_by(Metric) %>%
    kruskal_test(Value ~ treatment) %>%
    ungroup()
  pairwise_results <- suppressWarnings(
    alpha_long %>%
      group_by(Metric) %>%
      pairwise_wilcox_test(Value ~ treatment, p.adjust.method = "BH") %>%
      ungroup()
  )
  write_csv(kw_results, file.path(results_dir, "Alpha_Kruskal_Wallis.csv"))
  write_csv(pairwise_results, file.path(results_dir, "Alpha_Pairwise_Wilcoxon.csv"))

  # Faceted overview with the Kruskal-Wallis p-value per panel
  kw_labels <- kw_results %>%
    mutate(label = paste0("Kruskal-Wallis, p = ", signif(p, 2)))

  alpha_plot <- ggplot(alpha_long, aes(treatment, Value, fill = treatment)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15) +
    geom_text(data = kw_labels, aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.5, size = 3.2) +
    facet_wrap(~Metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = TREATMENT_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +  # room for the p-value
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(
      title    = paste(analysis_label, "- Functional Alpha Diversity"),
      subtitle = "MetaCyc pathways - PICRUSt2",
      x = "Treatment", y = NULL, fill = "Treatment"
    )
  save_plot(alpha_plot, results_dir, "Functional_Alpha_Diversity.png", 12, 5)

  # Shannon with pairwise brackets (significant pairs only)
  sig_pairs <- pairwise_results %>%
    filter(Metric == "Shannon", !is.na(p.adj), p.adj < ALPHA_CUTOFF)

  shannon_plot <- ggplot(alpha_df, aes(treatment, Shannon, fill = treatment)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15) +
    scale_fill_manual(values = TREATMENT_COLORS) +
    theme_pub() +
    labs(
      title    = paste(analysis_label, "- Functional Shannon Diversity"),
      subtitle = "PICRUSt2-predicted MetaCyc pathway profiles",
      x = "Treatment", y = "Shannon Diversity Index", fill = "Treatment"
    )

  if (nrow(sig_pairs) > 0) {
    y_max <- max(alpha_df$Shannon, na.rm = TRUE)
    step  <- diff(range(alpha_df$Shannon, na.rm = TRUE)) * 0.08
    sig_pairs <- sig_pairs %>% mutate(y.position = y_max + step * row_number())
    shannon_plot <- shannon_plot +
      stat_pvalue_manual(sig_pairs, label = "p.adj.signif", tip.length = 0.01) +
      labs(caption = "BH-adjusted pairwise Wilcoxon; only p ≤ 0.05 shown. * ≤ 0.05, ** ≤ 0.01, *** ≤ 0.001, **** ≤ 0.0001")
  } else {
    shannon_plot <- shannon_plot +
      labs(caption = "No pairwise differences at BH-adjusted p ≤ 0.05")
  }
  save_plot(shannon_plot, results_dir, "Shannon_Diversity.png", 8, 6)
}

# ---- Fold-change plots ----------------------------------------------------------

plot_fc_heatmap <- function(fc_data, results_dir, analysis_label) {
  top_features <- fc_data %>%
    group_by(feature) %>%
    summarise(max_abs_fc = max(abs(log2_fc)), .groups = "drop") %>%
    slice_max(max_abs_fc, n = TOP_N_HEATMAP, with_ties = FALSE)

  plot_data <- fc_data %>%
    filter(feature %in% top_features$feature) %>%
    complete(nesting(feature, description), comparison)  # grey = not significant

  p <- ggplot(plot_data, aes(comparison, description, fill = log2_fc)) +
    geom_tile(color = "grey85") +
    scale_fill_diverging(na.value = "grey90") +
    theme_pub() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text.y = element_text(size = 7)) +
    labs(
      title    = paste(analysis_label, "- Differentially Abundant Pathways"),
      subtitle = paste0("Top ", TOP_N_HEATMAP,
                        " pathways ranked by absolute Log2 fold change"),
      x = "Comparison", y = "MetaCyc Pathway", fill = "Log2 FC"
    )
  save_plot(p, results_dir, "fold_change_heatmap.png", 10, 8)
}

plot_fc_bars <- function(fc_data, results_dir, analysis_label) {
  plot_data <- fc_data %>%
    group_by(comparison) %>%
    slice_max(abs(log2_fc), n = TOP_N_BARS, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      description = str_wrap(description, width = 30),
      # reorder_within keeps bars sorted inside each facet even when the same
      # pathway appears in several comparisons
      description = reorder_within(description, log2_fc, comparison)
    )

  p <- ggplot(plot_data, aes(log2_fc, description, fill = log2_fc)) +
    geom_col() +
    facet_wrap(~comparison, scales = "free_y", ncol = 3) +
    scale_fill_diverging() +
    scale_y_reordered(expand = expansion(mult = c(0.03, 0.1))) +
    theme_pub() +
    theme(
      axis.text.y     = element_text(size = 7.5, lineheight = 0.85),
      strip.text      = element_text(face = "bold", size = 11,
                                     margin = margin(t = 5, b = 5)),
      panel.spacing.x = unit(1.4, "lines"),
      panel.spacing.y = unit(1.6, "lines"),
      plot.margin     = margin(10, 14, 10, 10)
    ) +
    labs(
      title = paste0(analysis_label, " - Top ", TOP_N_BARS,
                     " Differentially Abundant Pathways per Comparison"),
      x = "Log2 fold change", y = NULL, fill = "Log2 FC"
    )
  save_plot(p, results_dir, "Top15_Pathways_Per_Comparison.png", 22, 16)
}

# ---- Per-treatment z-score heatmap ----------------------------------------------

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
}

plot_zscore_heatmap <- function(daa_annotated, abundance_matrix, metadata,
                                results_dir, analysis_label) {
  long <- abundance_matrix %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "abundance") %>%
    left_join(metadata, by = "sample") %>%
    mutate(treatment = as.character(treatment))

  # Each treatment gets its OWN top pathways: lowest adjusted p among the
  # comparisons that treatment takes part in, and detected in that treatment.
  detected <- long %>%
    group_by(treatment, feature) %>%
    summarise(max_abundance = max(abundance, na.rm = TRUE), .groups = "drop") %>%
    filter(max_abundance > 0)

  top_features <- daa_annotated %>%
    filter(p_adjust < ALPHA_CUTOFF, !is.na(description)) %>%
    pivot_longer(c(group1, group2), names_to = "role", values_to = "treatment") %>%
    group_by(treatment, feature, description) %>%
    summarise(min_p = min(p_adjust), .groups = "drop") %>%
    semi_join(detected, by = c("treatment", "feature")) %>%
    group_by(treatment) %>%
    slice_min(min_p, n = TOP_N_ZSCORE, with_ties = FALSE) %>%
    ungroup()

  plot_data <- long %>%
    inner_join(top_features %>% select(treatment, feature, description),
               by = c("treatment", "feature")) %>%
    group_by(treatment, feature) %>%
    mutate(z_score = zscore(abundance)) %>%
    ungroup() %>%
    mutate(treatment = factor(treatment, levels = levels(metadata$treatment)))

  p <- ggplot(plot_data, aes(factor(day), description, fill = z_score)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_wrap(~treatment, ncol = 2, scales = "free") +
    scale_fill_diverging() +
    theme_pub() +
    theme(
      panel.grid    = element_blank(),
      axis.text.x   = element_text(angle = 0, size = 8),
      axis.text.y   = element_text(size = 7),
      panel.spacing = unit(1.2, "lines")
    ) +
    labs(
      title    = paste(analysis_label, "- Within-Treatment Patterns of Differentially Abundant Pathways"),
      subtitle = paste0("Top ", TOP_N_ZSCORE,
                        " significant MetaCyc pathways per treatment"),
      x = "Day", y = "MetaCyc Pathway", fill = "Z-score"
    )
  save_plot(p, results_dir, "zscore_heatmap.png", 18, 10)
}

# ------------------------------------------------------------------------------
# 5. Main workflow (same steps for both datasets)
# ------------------------------------------------------------------------------

run_picrust_analysis <- function(pathway_abundance, metadata, results_dir,
                                 analysis_label) {
  cat("\n---------------------------------\n", analysis_label,
      "\n---------------------------------\n", sep = "")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  abundance_matrix <- pathway_abundance %>%
    as.data.frame() %>%
    column_to_rownames(colnames(pathway_abundance)[1])

  # Differential abundance
  daa_annotated <- run_all_daa(abundance_matrix, metadata, comparisons)
  significant_pathways <- daa_annotated %>%
    filter(p_adjust < ALPHA_CUTOFF) %>%
    arrange(p_adjust)
  write_csv(daa_annotated,        file.path(results_dir, "all_daa_results.csv"))
  write_csv(significant_pathways, file.path(results_dir, "significant_pathways.csv"))
  cat(nrow(significant_pathways), "significant pathway hits\n")

  # Sample-level plots (do not depend on significant hits)
  plot_pca(abundance_matrix, metadata, results_dir, analysis_label)
  run_alpha_diversity(abundance_matrix, metadata, results_dir, analysis_label)

  # Plots that need significant hits
  fc_data <- prep_fold_changes(daa_annotated)
  if (nrow(fc_data) == 0) {
    message("No significant pathways for ", analysis_label,
            "; skipping fold-change and z-score plots.")
  } else {
    plot_fc_heatmap(fc_data, results_dir, analysis_label)
    plot_fc_bars(fc_data, results_dir, analysis_label)
    plot_zscore_heatmap(daa_annotated, abundance_matrix, metadata,
                        results_dir, analysis_label)
  }

  invisible(list(daa_annotated = daa_annotated,
                 significant_pathways = significant_pathways))
}

# ------------------------------------------------------------------------------
# 6. Run both analyses
# ------------------------------------------------------------------------------

results_all <- file.path(results_root, "All_MetaCyc")
results_org <- file.path(results_root, "Organonitrogen")
dir.create(results_org, recursive = TRUE, showWarnings = FALSE)

write_tsv(organonitrogen_pathways %>% select(feature, description),
          file.path(results_org, "organonitrogen_pathways.txt"))

all_results <- run_picrust_analysis(pathway_abundance, metadata,
                                    results_all, "All MetaCyc")

org_results <- run_picrust_analysis(organonitrogen_abundance, metadata,
                                    results_org, "Organonitrogen")
