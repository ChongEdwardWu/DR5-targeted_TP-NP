#!/usr/bin/env Rscript
# scRNA-seq — Step 2: Clustering, marker discovery, and annotation

# ------------------------- Setup & parameters ------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(future)
  library(stringr)
  library(SeuratObject)
  library(tibble)
  library(openxlsx2)
})

set.seed(123)
NWORKERS              <- as.integer(Sys.getenv("NWORKERS", unset = "20"))
PLAN_STRATEGY         <- Sys.getenv("PLAN", unset = "sequential")  # or "multicore"
PROJECT_BASE          <- Sys.getenv("PROJECT_BASE", unset = ".")
WORKDIR               <- Sys.getenv("WORKDIR", unset = file.path(PROJECT_BASE, "03_R"))
RES_DIR               <- file.path(WORKDIR, "results")
FIG_DIR               <- file.path(WORKDIR, "figures")
INTEGRATED_RDS        <- Sys.getenv("INTEGRATED_RDS", unset = file.path(RES_DIR, "01_Integration.rds"))
INCLUDE_DATE_IN_NAMES <- as.logical(Sys.getenv("INCLUDE_DATE_IN_NAMES", unset = "FALSE"))

if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
setwd(WORKDIR)

plan(PLAN_STRATEGY, workers = NWORKERS)

# ------------------------- Load integrated object -------------------------
seu <- readRDS(INTEGRATED_RDS)

# ------------------------- Section 1:  markers ---------------------
# Minimal lineage panel for quick sanity check via DotPlot.
marker_sets <- list(
  lineage_genes = c(
    # T cells
    "Cd3e",
    # CD8 T cells
    "Cd8a",
    # CD4 T cells
    "Cd4",
    # Treg
    "FoxP3",
    # B cell
    "Cd79a",
    # GC B cells
    "Fas",
    # NK
    "Ncr1", "Klrb1c",
    # DC
    "Ftl3", "Itgax"
  )
)

gene_panel <- unique(unlist(marker_sets_l1))
gene_panel <- gene_panel[gene_panel %in% rownames(seu)]

try({
  p <- VlnPlot(seu, assay = "SCT", features = gene_panel, stack = T, pt.size=0.1)
  ggsave(file.path(FIG_DIR, "Step2_vlnplot.png"), p, width = 180, height = 160, units = "mm", dpi = 150, bg = "white")
}, silent = TRUE)

# Coarse annotation (cluster-level). Keep the original mapping.
seu$CellType <- factor(dplyr::recode(
  seu$seurat_clusters,
  "1"="Bcell","2"="CD8T","3"="CD4T","4"="CD8T","5"="Treg","6"="Bcell",
  "7"="CD8T","8"="NK","9"="DC", "10"="GC_B"
))

# ------------------------- Save annotated object & session info ------------
res_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("02_Annotation%s.rds", format(Sys.Date(), "%Y%m%d")) else "02_Annotation.rds"
saveRDS(seu, file = file.path(RES_DIR, res_file))

sess_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("session_info_step2_%s.txt", format(Sys.Date(), "%Y%m%d")) else "session_info_step2.txt"
writeLines(capture.output(sessionInfo()), file.path(RES_DIR, sess_file))
