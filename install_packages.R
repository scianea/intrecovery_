# Installs every package used by the PICRUSt2 and RDA scripts
cran_pkgs <- c("tidyverse", "ggpicrust2", "vegan", "rstatix", "ggpubr", "tidytext",
               "ragg", "patchwork", "readxl", "ggrepel")
bioc_pkgs <- c("DESeq2")

missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing_cran) > 0) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE)

message("All packages installed.")
