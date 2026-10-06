# ==============================================================================
# v2 - same figure as z_score.R, built from the v2 DESeq2 results
# (picrust2_dual_analysis_v2). The original script and figures are kept unchanged.
# Run picrust2_analysis/01_picrust2_dual_analysis.R first.
# v3 (cosmetic): one-line pathway labels so they no longer overlap, HTML codes in
# names turned into symbols (&beta; -> β), a few long names abbreviated (SP, GlcNAc ...).
# ==============================================================================
# COMBINED Z-SCORE FIGURE (journal-ready)
#   A = All MetaCyc pathways
#   B = Organonitrogen pathways
# Reads the results already written by 01_picrust2_dual_analysis.R
# (all_daa_results.csv in each folder), so DESeq2 does NOT need to be re-run.
# No titles; bold "A"/"B" in the top-left corner of each panel; one shared
# colour bar; sizes are set in points for the final printed width.
# ==============================================================================

library(tidyverse)
library(patchwork)
library(ragg)

# ------------------------------------------------------------------------------
# 1. Config
# ------------------------------------------------------------------------------

picrust2_dir <- "data/picrust2_output"
results_root <- "results/picrust2_analysis"   # output of 01_picrust2_dual_analysis.R
out_dir      <- file.path(results_root, "Combined_Figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

daa_files <- c(
  A = file.path(results_root, "All_MetaCyc",    "all_daa_results.csv"),
  B = file.path(results_root, "Organonitrogen", "all_daa_results.csv")
)

ALPHA_CUTOFF     <- 0.05
TOP_N_ZSCORE     <- 10    # pathways per treatment (the standalone figures used 15)
LABEL_SIZE       <- 6     # pt; pathway names on the y axes (one line each, never wrapped)
LABEL_MODE       <- "name"  # "name" = pathway names (shortened below); "id" = MetaCyc IDs
# (narrow labels, wider heat maps; full names stay in the data CSV)

# Final printed size in inches. 7.5 in = full page width for most journals;
# check your journal's limit and edit these two numbers if needed.
FIG_WIDTH  <- 7.5
FIG_HEIGHT <- 9.5

# ------------------------------------------------------------------------------
# 2. Theme (same theme_pub as every other manuscript figure, print-sized)
# ------------------------------------------------------------------------------

PUB_FONT      <- "Times New Roman"
PUB_BASE_SIZE <- 9        # pt; all other sizes below scale from this

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

# ------------------------------------------------------------------------------
# 3. Data: abundance + metadata (same parsing as the main PICRUSt2 script)
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

abundance_matrix <- pathway_abundance %>%
  as.data.frame() %>%
  column_to_rownames(colnames(pathway_abundance)[1])

metadata <- parse_sample_name(colnames(abundance_matrix)) %>%
  mutate(treatment = factor(treatment, levels = paste0("R", 1:6)))

bad_samples <- metadata %>% filter(is.na(treatment) | is.na(day)) %>% pull(sample)
if (length(bad_samples) > 0) {
  stop("Samples not matched by parse_sample_name(): ", paste(bad_samples, collapse = ", "))
}

# ------------------------------------------------------------------------------
# 4. Z-score data (same selection logic as plot_zscore_heatmap() in the main script)
# ------------------------------------------------------------------------------

# One-line axis labels. The full names are kept in the `description` column of the data
# CSV. Abbreviations used here must be defined in the figure legend:
#   SP, superpathway of; GlcNAc, N-acetylglucosamine; ManNAc, N-acetylmannosamine;
#   Neu5Ac, N-acetylneuraminate.
# To keep a name in full, delete its line from the replacement list below.
GREEK <- c("&alpha;" = "α", "&beta;" = "β", "&gamma;" = "γ",
           "&delta;" = "δ", "&epsilon;" = "ε", "&omega;" = "ω")

shorten_label <- function(x) {
  x <- str_replace_all(x, GREEK)
  str_replace_all(x, c(
    "superpathay"              = "superpathway",   # typo in the MetaCyc description
    "^superpathway of "        = "SP: ",
    "N-acetylglucosamine"      = "GlcNAc",
    "N-acetylmannosamine"      = "ManNAc",
    "N-acetylneuraminate"      = "Neu5Ac",
    "Saccharomyces cerevisiae" = "S. cerevisiae"
  ))
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
}

build_zscore_data <- function(daa_annotated) {
  long <- abundance_matrix %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "abundance") %>%
    left_join(metadata, by = "sample") %>%
    mutate(treatment = as.character(treatment))
  
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
  
  long %>%
    inner_join(top_features %>% select(treatment, feature, description),
               by = c("treatment", "feature")) %>%
    group_by(treatment, feature) %>%
    mutate(z_score = zscore(abundance)) %>%
    ungroup() %>%
    mutate(
      treatment = factor(treatment, levels = levels(metadata$treatment)),
      label     = if (LABEL_MODE == "id") feature else shorten_label(description)
    )
}

zscore_data <- map(daa_files, ~ build_zscore_data(read_csv(.x, show_col_types = FALSE)))

# One symmetric colour range for BOTH panels so a single colour bar is valid
z_limit <- max(abs(unlist(map(zscore_data, "z_score"))), na.rm = TRUE)

# Supplementary table: full pathway names (description) and the axis label used (label)
zscore_data %>%
  imap(~ mutate(.x, panel = .y)) %>%
  bind_rows() %>%
  write_csv(file.path(out_dir, "Figure_ZScore_Combined_data.csv"))

# ------------------------------------------------------------------------------
# 5. One panel (no title, no subtitle)
# ------------------------------------------------------------------------------

plot_zscore_panel <- function(plot_data) {
  ggplot(plot_data, aes(factor(day), label, fill = z_score)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_wrap(~treatment, ncol = 2, scales = "free") +
    scale_fill_gradient2(low = "red", mid = "white", high = "darkblue",
                         midpoint = 0, limits = c(-z_limit, z_limit)) +
    labs(x = "Day", y = if (LABEL_MODE == "id") "MetaCyc pathway ID" else "MetaCyc pathway",
         fill = "Z-score") +
    theme_pub() +
    theme(
      panel.grid    = element_blank(),
      axis.text.x   = element_text(size = 6.5),
      axis.text.y   = element_text(size = LABEL_SIZE),
      panel.spacing = unit(0.6, "lines"),
      plot.margin   = margin(t = 16, r = 6, b = 4, l = 16)  # room for the A/B tag
    )
}

panel_A <- plot_zscore_panel(zscore_data[["A"]])
panel_B <- plot_zscore_panel(zscore_data[["B"]])

# ------------------------------------------------------------------------------
# 6. Combine: A over B, shared colour bar at the bottom, bold A/B top-left
# ------------------------------------------------------------------------------

combined <- (panel_A / panel_B) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position     = "bottom",
    plot.tag            = element_text(family = PUB_FONT, face = "bold", size = 16),
    plot.tag.position   = "topleft"
  )

# ------------------------------------------------------------------------------
# 7. Save: PNG + TIFF (600 dpi) and a vector PDF
# ------------------------------------------------------------------------------

base_name <- file.path(out_dir, "Figure_ZScore_Combined")

ggsave(paste0(base_name, ".png"),  combined, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = 600)
ggsave(paste0(base_name, ".tiff"), combined, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = 600,
       compression = "lzw")

tryCatch(
  ggsave(paste0(base_name, ".pdf"), combined, width = FIG_WIDTH, height = FIG_HEIGHT,
         device = cairo_pdf),
  error = function(e) message("PDF not written: ", conditionMessage(e))
)

message("Saved to: ", out_dir)