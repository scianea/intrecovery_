# ============================================================
# RDA Analysis - Full Pipeline (v4, 3009)
# Runs the SAME RDA on BOTH PICRUSt2 datasets so they can be compared:
#   1. All MetaCyc pathways      (PC1/PC2 from results/picrust2_analysis/All_MetaCyc)
#   2. Organonitrogen pathways   (PC1/PC2 from results/picrust2_analysis/Organonitrogen)
# Explanatory variable: Treatment (R1-R6)
# Response variables:   chemistry + community (EMU PCoA1/2) + function (PICRUSt2 PC1/2)
# Only PC1/PC2 differ between the two analyses.
# ============================================================

library(readxl); library(dplyr); library(tidyr); library(stringr)
library(readr); library(purrr); library(tibble); library(vegan)
library(ggplot2); library(ggrepel); library(patchwork)

# ------------------------------------------------------------------------------
# 0. Publication font & theme (Times New Roman, one consistent size scale)
# ------------------------------------------------------------------------------
PUB_FONT <- "Times New Roman"
PUB_BASE_SIZE <- 12

theme_pub <- function(base_size = PUB_BASE_SIZE, base_family = PUB_FONT) {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title = element_text(
        family = base_family, face = "bold",
        size = base_size + 3, hjust = 0.5,
        margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        family = base_family, size = base_size + 1, hjust = 0.5,
        margin = margin(b = 6)
      ),
      plot.caption = element_text(family = base_family, size = base_size - 3),
      axis.title = element_text(family = base_family, size = base_size),
      axis.text = element_text(family = base_family, size = base_size - 2),
      legend.title = element_text(family = base_family, size = base_size),
      legend.text = element_text(family = base_family, size = base_size - 1),
      strip.text = element_text(family = base_family, size = base_size, face = "bold")
    )
}
theme_set(theme_pub())

# Locked treatment colours: identical to TREATMENT_COLORS in
# emu_final_analysis_3009.R and picrust2_first_manuscript_edited_3009.R
TREATMENT_COLORS <- c(
  "R1" = "#F8766D",  # Control
  "R2" = "#B79F00",  # BES
  "R3" = "#00BA38",  # Cysteine
  "R4" = "#00BFC4",  # Ascorbate
  "R5" = "#619CFF",  # Titanium
  "R6" = "#F564E3"   # Mix
)

# ------------------------------------------------------------------------------
# 1. Paths (relative to the repository root)
# ------------------------------------------------------------------------------
input_dir       <- "data/chemistry"                      # chemistry workbook (All_Results_RDA_M1.xlsx)
output_root     <- "results/rda_analysis"                # one subfolder per analysis
base_dir        <- "data/emu_results"                    # raw EMU rel-abundance files
emu_results_dir <- file.path(base_dir, "final_analysis")  # selected_best_runs.csv + PCoA_coordinates.csv

picrust_dirs <- c(
  "All MetaCyc"    = "results/picrust2_analysis/All_MetaCyc",
  "Organonitrogen" = "results/picrust2_analysis/Organonitrogen"
)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
set.seed(123)  # permutation tests (anova.cca) are reproducible and comparable

# ------------------------------------------------------------
# MASTER TREATMENT RECODER - the ONE place treatment names get standardized
# Accepts long names (Control, BES, ...) AND the R1-R6 codes now written by
# the EMU and PICRUSt2 scripts, so all sample keys use the same naming.
# ------------------------------------------------------------
standardize_treatment_col <- function(x) {
  x <- str_squish(str_replace_all(x, "\u00A0", " "))  # kill non-breaking spaces
  x <- dplyr::recode(x, R1 = "Control", R2 = "BES", R3 = "Cysteine",
                     R4 = "Ascorbate", R5 = "Titanium", R6 = "Mix")
  dplyr::case_when(
    str_detect(x, regex("^Ascorbat|^Ascorbic", ignore_case = TRUE)) ~ "Ascorbate",
    str_detect(x, regex("^Control.*BES|^Bes",   ignore_case = TRUE)) ~ "BES",
    str_detect(x, regex("^Control$",            ignore_case = TRUE)) ~ "Control",
    str_detect(x, regex("^Cysteine",            ignore_case = TRUE)) ~ "Cysteine",
    str_detect(x, regex("^TiCl|^Titanium",      ignore_case = TRUE)) ~ "Titanium",
    str_detect(x, regex("^Mixture|^Mixed|^Mix", ignore_case = TRUE)) ~ "Mix",
    TRUE ~ x  # unmatched stays visible, never silently becomes NA
  )
}

make_sample_key <- function(treatment, day) paste(treatment, day, sep = "_")

# ------------------------------------------------------------
# PART 1: Environmental / chemistry data
# ------------------------------------------------------------
file_path <- file.path(input_dir, "All_Results_RDA_M1.xlsx")

nh4_nh3   <- read_excel(file_path, sheet = "NH4_NH3")
ph_orp    <- read_excel(file_path, sheet = "PH_ORP")
aa_amines <- read_excel(file_path, sheet = "AA_Amines")
scfas     <- read_excel(file_path, sheet = "SCFAs_Total C")
names(ph_orp) <- str_trim(names(ph_orp))

nh4_nh3   <- nh4_nh3   %>% dplyr::mutate(Treatment = standardize_treatment_col(Treatment))
ph_orp    <- ph_orp    %>% dplyr::mutate(Treatment = standardize_treatment_col(Treatment))
aa_amines <- aa_amines %>% dplyr::mutate(Treatment = standardize_treatment_col(Treatment))
scfas     <- scfas     %>% dplyr::mutate(Treatment = standardize_treatment_col(Treatment))

ph_orp <- ph_orp %>%
  dplyr::mutate(pH = as.numeric(str_replace(pH, ",", ".")),
                `ORP (mV)` = as.numeric(str_replace(`ORP (mV)`, ",", ".")))

nh4_nh3 <- nh4_nh3 %>% dplyr::mutate(Sample = if_else(Sample == "CY1R2", "CY2", Sample))

nh4_nh3   <- nh4_nh3   %>% dplyr::rename(RetentionTime = `Retention Time (day)`)
ph_orp    <- ph_orp    %>% dplyr::rename(RetentionTime = `Retention Time (day)`)
aa_amines <- aa_amines %>% dplyr::rename(RetentionTime = Day)
scfas     <- scfas     %>% dplyr::rename(RetentionTime = Day)

nh4_nh3   <- nh4_nh3   %>% dplyr::mutate(sample_key = make_sample_key(Treatment, RetentionTime))
ph_orp    <- ph_orp    %>% dplyr::mutate(sample_key = make_sample_key(Treatment, RetentionTime))
aa_amines <- aa_amines %>% dplyr::mutate(sample_key = make_sample_key(Treatment, RetentionTime))
scfas     <- scfas     %>% dplyr::mutate(sample_key = make_sample_key(Treatment, RetentionTime))

env_merged <- nh4_nh3 %>%
  dplyr::select(sample_key, Treatment, RetentionTime,
                NH4_mean = `NH4_mean (mmol/L)`, NH3_mean = `NH3_mean (mmol/L)`) %>%
  dplyr::left_join(aa_amines %>% dplyr::select(sample_key,
                                               TotalAmine = `Total Amine Production (mM)`,
                                               TotalAminoAcid = `Total Amino Acid Production (mM)`,
                                               TotalOrganoN = `Total Organonitrogen (mM)`), by = "sample_key") %>%
  dplyr::left_join(scfas %>% dplyr::select(sample_key, TotalSCFA = `Total SCFA Concentration (mM)`),
                   by = "sample_key") %>%
  dplyr::left_join(ph_orp %>% dplyr::select(sample_key, pH, `ORP (mV)`), by = "sample_key")

dim(env_merged)              # expect 59 x 11
colSums(is.na(env_merged))   # expect NH4_mean/NH3_mean = 4 each, rest 0

# ------------------------------------------------------------
# PART 2: Community data (EMU species counts)
# Kept here for provenance/future use, even though the final model below
# uses PCoA1/PCoA2 (already-reduced compositional summary) rather than the
# full species matrix as its response.
# ------------------------------------------------------------
run_to_folder <- c("Batch1Final" = "Batch1_Final", "TiCtMix" = "TiCtMix",
                   "CysAsc" = "", "CtrlBES" = "")

selected_runs <- read_csv(file.path(emu_results_dir, "selected_best_runs.csv"),
                          show_col_types = FALSE) %>%
  dplyr::mutate(
    Treatment = standardize_treatment_col(treatment),
    sample_key = make_sample_key(Treatment, day),
    folder = run_to_folder[run],
    file_name = paste0(run, "_", barcode, "_merged.fastq_rel-abundance.tsv"),
    file_path = file.path(base_dir, folder, file_name)
  )

read_one_emu <- function(path, key) {
  read_tsv(path, show_col_types = FALSE) %>%
    dplyr::select(species, `estimated counts`) %>%
    dplyr::mutate(sample_key = key)
}
all_emu <- map2_dfr(selected_runs$file_path, selected_runs$sample_key, read_one_emu)

community_matrix <- all_emu %>%
  dplyr::filter(!is.na(species)) %>%
  dplyr::group_by(sample_key, species) %>%
  dplyr::summarise(count = sum(`estimated counts`, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = species, values_from = count, values_fill = 0) %>%
  tibble::column_to_rownames("sample_key") %>%
  as.matrix()

dim(community_matrix)  # expect 59 x ~1300

# ------------------------------------------------------------
# PART 3: EMU PCoA1/PCoA2 (shared by both analyses) + derived variables
# ------------------------------------------------------------
pcoa <- read_csv(file.path(emu_results_dir, "PCoA_coordinates.csv"),
                 show_col_types = FALSE) %>%
  dplyr::mutate(Treatment = standardize_treatment_col(treatment),
                sample_key = make_sample_key(Treatment, day))

env_merged <- env_merged %>%
  dplyr::left_join(pcoa %>% dplyr::select(sample_key, PCoA1, PCoA2), by = "sample_key")

# --- Derive TAN and FAN ---
env_merged <- env_merged %>%
  dplyr::mutate(
    TAN = NH4_mean + NH3_mean,
    FAN = NH3_mean
  )

# --- Final alignment check (chemistry vs community rows) ---
env_merged <- env_merged %>% dplyr::arrange(sample_key)
community_matrix <- community_matrix[order(rownames(community_matrix)), ]

identical(env_merged$sample_key, rownames(community_matrix))  # must be TRUE
sum(is.na(env_merged$PCoA1))  # must be 0
sum(is.na(env_merged$PCoA2))  # must be 0

# --- Relabel treatments: R1-R6 ---
relabel_treatment <- function(x) {
  dplyr::recode(as.character(x),
                "Control"   = "R1",
                "BES"       = "R2",
                "Cysteine"  = "R3",
                "Ascorbate" = "R4",
                "Titanium"  = "R5",
                "Mix"       = "R6"
  )
}

env_merged <- env_merged %>%
  dplyr::mutate(Treatment = factor(relabel_treatment(Treatment),
                                   levels = c("R1", "R2", "R3", "R4", "R5", "R6")))

# ------------------------------------------------------------
# PART 4: Triplot helper (used for RDA1 vs RDA2 and RDA1 vs RDA3)
# Styled like the PCoA/PCA: filled 95% ellipses (alpha 0.2), locked colours.
# ------------------------------------------------------------
format_p <- function(p) {
  if (is.na(p)) "p = NA" else if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
}

make_triplot <- function(rda_model, env_complete, x_axis, y_axis,
                         anova_overall, anova_axis) {
  ax      <- paste0("RDA", c(x_axis, y_axis))
  choices <- c(x_axis, y_axis)
  imp     <- summary(rda_model)$cont$importance
  
  site_scores <- as.data.frame(scores(rda_model, display = "sites", choices = choices))
  site_scores$sample_key <- env_complete$sample_key
  site_scores <- site_scores %>%
    dplyr::left_join(env_complete %>% dplyr::select(sample_key, Treatment), by = "sample_key")
  
  response_scores <- as.data.frame(scores(rda_model, display = "species", choices = choices))
  response_scores$label <- rownames(response_scores)
  
  centroid_scores <- as.data.frame(scores(rda_model, display = "cn", choices = choices))
  centroid_scores$label <- stringr::str_remove(rownames(centroid_scores), "^Treatment")
  
  # Permutation-test results shown on the figure
  p_axis_x <- anova_axis$`Pr(>F)`[x_axis]
  p_axis_y <- anova_axis$`Pr(>F)`[y_axis]
  test_caption <- sprintf(
    "Permutation test (999 permutations): pseudo-F = %.2f, %s; adjusted R\u00B2 = %.2f",
    anova_overall$F[1], format_p(anova_overall$`Pr(>F)`[1]),
    RsquareAdj(rda_model)$adj.r.squared
  )
  
  arrow_mult <- 0.8 * max(abs(as.matrix(site_scores[, ax]))) /
    max(abs(as.matrix(response_scores[, ax])))
  
  ggplot() +
    geom_point(data = site_scores,
               aes(x = .data[[ax[1]]], y = .data[[ax[2]]], color = Treatment),
               size = 3) +
    stat_ellipse(data = site_scores,
                 aes(x = .data[[ax[1]]], y = .data[[ax[2]]], fill = Treatment),
                 geom = "polygon", alpha = 0.2, color = NA) +
    geom_point(data = centroid_scores,
               aes(x = .data[[ax[1]]], y = .data[[ax[2]]]),
               shape = 4, size = 2, stroke = 2, color = "black") +
    geom_text_repel(data = centroid_scores,
                    aes(x = .data[[ax[1]]], y = .data[[ax[2]]], label = label),
                    color = "black", fontface = "bold", size = 4, family = PUB_FONT,
                    box.padding = 1.2, point.padding = 0.5,
                    force = 8, max.overlaps = Inf, min.segment.length = 0) +
    geom_segment(data = response_scores,
                 aes(x = 0, y = 0,
                     xend = .data[[ax[1]]] * arrow_mult,
                     yend = .data[[ax[2]]] * arrow_mult),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "darkred", linewidth = 0.6) +
    geom_text_repel(data = response_scores,
                    aes(x = .data[[ax[1]]] * arrow_mult,
                        y = .data[[ax[2]]] * arrow_mult, label = label),
                    color = "darkred", fontface = "italic", size = 3.5, family = PUB_FONT) +
    scale_color_manual(values = TREATMENT_COLORS) +
    scale_fill_manual(values = TREATMENT_COLORS) +
    labs(
      x = paste0(ax[1], " (", round(imp[2, x_axis] * 100, 1), "%, ", format_p(p_axis_x), ")"),
      y = paste0(ax[2], " (", round(imp[2, y_axis] * 100, 1), "%, ", format_p(p_axis_y), ")"),
      caption = test_caption,
      title = paste0("RDA", x_axis, " vs RDA", y_axis),
      color = "Treatment", fill = "Treatment"
    ) +
    theme_pub() +
    theme(
      legend.position = "right",
      plot.caption    = element_text(margin = margin(t = 14)),  # pushes the caption clear of the x-axis title
      plot.margin     = margin(t = 20, r = 8, b = 12, l = 20)  # top/left room for the A/B tags
    )
}

# ------------------------------------------------------------
# PART 5: One full RDA for one PICRUSt2 dataset
# ------------------------------------------------------------
run_rda_analysis <- function(analysis_label, picrust_dir, env) {
  out_dir <- file.path(output_root, str_replace_all(analysis_label, " ", "_"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\n=== RDA:", analysis_label, "===\n")
  
  # PICRUSt2 PC1/PC2 for this dataset
  picrust_pca <- read_csv(file.path(picrust_dir, "PCA_PC1_PC2.csv"),
                          show_col_types = FALSE) %>%
    dplyr::mutate(Treatment = standardize_treatment_col(treatment),
                  sample_key = make_sample_key(Treatment, day))
  
  env <- env %>%
    dplyr::left_join(picrust_pca %>% dplyr::select(sample_key, PC1, PC2), by = "sample_key")
  cat("Samples without PICRUSt2 PC scores:", sum(is.na(env$PC1)), "(expect 0)\n")
  
  env_response <- env %>%
    dplyr::select(pH, `ORP (mV)`, TotalSCFA, TAN, FAN,
                  TotalAminoAcid, TotalAmine, PCoA1, PCoA2, PC1, PC2) %>%
    as.data.frame()
  rownames(env_response) <- env$sample_key
  
  # Complete cases only (no NAs across the response variables)
  complete_idx <- complete.cases(env_response)
  cat("Rows dropped for missing values:", sum(!complete_idx), "(expect 4)\n")
  env_response_complete <- env_response[complete_idx, ]
  env_complete          <- env[complete_idx, ]
  
  # Treatment as explanatory variable, everything else as response
  rda_model <- rda(env_response_complete ~ Treatment, data = env_complete, scale = TRUE)
  
  anova_overall <- anova.cca(rda_model, permutations = 999)                # whole model
  anova_terms   <- anova.cca(rda_model, by = "terms", permutations = 999)  # Treatment
  anova_axis    <- anova.cca(rda_model, by = "axis", permutations = 999)   # per RDA axis
  
  capture.output(
    cat("RDA:", analysis_label, "\n\n"),
    print(summary(rda_model)),
    cat("\n--- Overall model ---\n"), print(anova_overall),
    cat("\n--- By term ---\n"),       print(anova_terms),
    cat("\n--- By axis ---\n"),       print(anova_axis),
    file = file.path(out_dir, "RDA_summary_and_tests.txt")
  )
  write_csv(env_complete, file.path(out_dir, "RDA_input_data.csv"))
  
  # Triplots
  p12 <- make_triplot(rda_model, env_complete, 1, 2, anova_overall, anova_axis)
  p13 <- make_triplot(rda_model, env_complete, 1, 3, anova_overall, anova_axis)
  
  ggsave(file.path(out_dir, "rda1v2.png"), p12, width = 10, height = 8, dpi = 300)
  ggsave(file.path(out_dir, "rda1v3.png"), p13, width = 10, height = 8, dpi = 300)
  
  # Stacked: RDA1 vs RDA2 on top, RDA1 vs RDA3 below, one shared legend
  # Bold A / B tags in the top-left corner of each panel
  combined_plot <- (p12 / p13) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(
      legend.position   = "right",
      plot.tag          = element_text(family = PUB_FONT, face = "bold", size = 20),
      plot.tag.position = "topleft"
    )
  ggsave(file.path(out_dir, "rda_triplot_combined.png"),
         combined_plot, width = 10, height = 16, dpi = 300)
  
  invisible(list(model = rda_model, anova_overall = anova_overall,
                 anova_terms = anova_terms, anova_axis = anova_axis,
                 n = nrow(env_response_complete)))
}

# ------------------------------------------------------------
# PART 6: Run both analyses and compare them side by side
# ------------------------------------------------------------
rda_results <- imap(picrust_dirs, ~ run_rda_analysis(.y, .x, env_merged))

rda_comparison <- imap_dfr(rda_results, function(r, label) {
  tibble(
    analysis            = label,
    n_samples           = r$n,
    constrained_percent = 100 * r$model$CCA$tot.chi / r$model$tot.chi,
    adjusted_R2         = RsquareAdj(r$model)$adj.r.squared,
    F_overall           = r$anova_overall$F[1],
    p_overall           = r$anova_overall$`Pr(>F)`[1]
  )
})

print(rda_comparison)
write_csv(rda_comparison, file.path(output_root, "RDA_comparison_All_vs_Organonitrogen.csv"))