#!/usr/bin/env Rscript
# scRNA-seq — Step 1: Integration & initial clustering (Seurat SCT)

# ------------------------- Setup & parameters ------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(patchwork)
  library(tidyverse)
  library(ggplot2)
  library(clustree)
  library(future)
  library(openxlsx2)
})

set.seed(123)
NWORKERS              <- as.integer(Sys.getenv("NWORKERS", unset = "8"))
PLAN_STRATEGY         <- Sys.getenv("PLAN", unset = "sequential")  # or "multicore"
PROJECT_BASE          <- Sys.getenv("PROJECT_BASE", unset = ".")
WORKDIR               <- Sys.getenv("WORKDIR", unset = file.path(PROJECT_BASE, "03_R"))
QC_DIR                <- Sys.getenv("QC_DIR", unset = file.path(WORKDIR, "QC"))
SAMPLES_ENV           <- Sys.getenv("SAMPLES", unset = "PBS,MD5,TP")
INCLUDE_DATE_IN_NAMES <- as.logical(Sys.getenv("INCLUDE_DATE_IN_NAMES", unset = "FALSE"))
RESULTS_BASENAME      <- Sys.getenv("RESULTS_BASENAME", unset = "01_Integration")
FIG_DIR               <- file.path(WORKDIR, "figures")
RES_DIR               <- file.path(WORKDIR, "results")

if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)
setwd(WORKDIR)

# Parse sample list (defaults to two example samples). You may also list.dirs(QC_DIR).
samples <- trimws(strsplit(SAMPLES_ENV, ",")[[1]])

# -------------------- Load per-sample Seurat objects -----------------------
Seu_list <- list()
for (s in samples) {
  rds_path <- file.path(QC_DIR, s, sprintf("%s_seu.rds", s))
  if (!file.exists(rds_path)) stop("Missing input: ", rds_path)
  Seu_list[[s]] <- readRDS(rds_path)
}

# -------------------- SCT Integration (RPCA) -------------------------------
features <- SelectIntegrationFeatures(object.list = Seu_list, nfeatures = 2500)
Seu_list <- PrepSCTIntegration(object.list = Seu_list, anchor.features = features, verbose = TRUE)

# Optional reference: set SEURAT_REF to a sample name to anchor to it
SEURAT_REF <- Sys.getenv("SEURAT_REF", unset = "MD5")
reference_idx <- if (nzchar(SEURAT_REF)) which(names(Seu_list) == SEURAT_REF) else NULL
anchors <- FindIntegrationAnchors(object.list = Seu_list, normalization.method = "SCT",
                                  reduction = "rpca", anchor.features = features,
                                  reference = reference_idx, verbose = TRUE)
seu <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE)

# -------------------- PCA/UMAP on curated features -------------------------
DefaultAssay(seu) <- "integrated"

# Build a list of features to exclude from PCA (histone/mito/ribo/etc.)
# Use RNA assay rownames if present, else fallback to integrated feature names.
rn <- if ("RNA" %in% names(seu@assays)) rownames(seu@assays$RNA@counts) else rownames(seu)
cc_genes      <- unique(unlist(Seurat::cc.genes))
hist_genes    <- grep("^Hist", rn, ignore.case = TRUE, value = TRUE)
hb_genes      <- grep("^Hb[ab]-|^HB[^(P)]", rn, value = TRUE)
mt_genes      <- grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l", rn, ignore.case = TRUE, value = TRUE)
rps_genes     <- grep("^Rp[sl]|^RP[SL]", rn, ignore.case = TRUE, value = TRUE)
rik_genes     <- grep("^Rik", rn, ignore.case = TRUE, value = TRUE)
alu_genes     <- grep("^AL", rn, ignore.case = TRUE, value = TRUE)
pseudo_genes  <- grep("-rs|-ps", rn, ignore.case = TRUE, value = TRUE)
mir_genes     <- grep("^Mir", rn, ignore.case = TRUE, value = TRUE)
gencode_genes <- grep("^Gm", rn, ignore.case = TRUE, value = TRUE)

bad_features <- unique(c(cc_genes, hist_genes, hb_genes, mt_genes, rps_genes,
                         rik_genes, alu_genes, pseudo_genes, mir_genes, gencode_genes))

PCA_features <- setdiff(seu@assays$integrated@var.features, bad_features)
seu <- RunPCA(seu, npcs = 50, verbose = TRUE, features = PCA_features)

# Elbow (save to file for reproducibility)
elbow_plot <- ElbowPlot(seu, ndims = 50)
ggsave(file.path(FIG_DIR, "Step1_elbow_plot.png"), elbow_plot, width = 140, height = 120, units = "mm", dpi = 150)

# UMAP/Neighbors (dims 1:15 as in your script)
plan(PLAN_STRATEGY, workers = NWORKERS)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:15)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:15)

# Quick UMAP by sample/group (if `group` meta present from previous steps)
umap_by_group <- try({
  DimPlot(seu, reduction = "umap", group.by = c("group"), label.size = 2)
}, silent = TRUE)
if (!inherits(umap_by_group, "try-error")) {
  ggsave(file.path(FIG_DIR, "Step1_umap_by_group.png"), umap_by_group, width = 160, height = 140, units = "mm", dpi = 150, bg = "white")
}

# -------------------- Clustering across resolutions & Filter bad cells/clusters -----------------------
seu <- FindClusters(seu, resolution = c(seq(0.1, 1.0, 0.1)))
# Relabel clusters from 1..K (avoid "0") per original practice
meta_cols <- grep("integrated_snn_res.", colnames(seu@meta.data), value = TRUE)
for (cn in meta_cols) {
  k <- seu@meta.data[[cn]]
  levels(k) <- as.character(seq_len(nlevels(k)))
  seu@meta.data[[cn]] <- k
}

# Clustree visualization across resolutions
ct <- clustree(seu@meta.data, prefix = "integrated_snn_res.", return = "plot")
ggsave(file.path(FIG_DIR, "Step1_clustree.png"), ct, width = 180, height = 140, units = "mm", dpi = 150)

# Choose a working resolution for downstream visualization/exports
RES_CHOSEN <- as.numeric(Sys.getenv("RES_CHOSEN", unset = "0.7"))
Idents(seu) <- seu$seurat_clusters <- factor(seu@meta.data[[paste0("integrated_snn_res.", RES_CHOSEN)]])

# UMAP split by sample, colored by chosen clusters
umap_split <- try({
  DimPlot(seu, reduction = "umap", split.by = "group",
          group.by = paste0("integrated_snn_res.", RES_CHOSEN), label = TRUE) +
    coord_fixed(ratio = 1)
}, silent = TRUE)
if (!inherits(umap_split, "try-error")) {
  ggsave(file.path(FIG_DIR, sprintf("Step1_umap_split_res%s.png", RES_CHOSEN)), umap_split,
         width = 180, height = 140, units = "mm", dpi = 150, bg = "white")
}

# QC violin plots (if the relevant columns exist)
qc_feats <- c("sum", "detected", "subsets_Mito_percent", "subsets_Rp_percent", "subsets_Heatshock_percent")
qc_feats <- qc_feats[qc_feats %in% colnames(seu@meta.data)]
if (length(qc_feats)) {
  vln <- VlnPlot(seu, features = qc_feats, pt.size = 0)
  ggsave(file.path(FIG_DIR, "Step1_qc_vln.png"), vln, width = 200, height = 140, units = "mm", dpi = 150, bg = "white")
}

# Filter bad cells/clusters
discardCl <- seu$seurat_clusters %in% c(8, 11)
seu_filt <- seu[, !(discardCl)]
seu <- seu_filt

seu <- RunPCA(seu, npcs = 50, verbose = TRUE, features = PCA_features)

# Elbow (save to file for reproducibility)
elbow_plot <- ElbowPlot(seu, ndims = 50)
ggsave(file.path(FIG_DIR, "Step1_elbow_plot.png"), elbow_plot, width = 140, height = 120, units = "mm", dpi = 150)

# UMAP/Neighbors (dims 1:20 as in your script)
plan(PLAN_STRATEGY, workers = NWORKERS)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:20)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)

seu <- FindClusters(seu, resolution = c(seq(0.1, 1.0, 0.1)))
# Relabel clusters from 1..K (avoid "0") per original practice
meta_cols <- grep("integrated_snn_res.", colnames(seu@meta.data), value = TRUE)
for (cn in meta_cols) {
  k <- seu@meta.data[[cn]]
  levels(k) <- as.character(seq_len(nlevels(k)))
  seu@meta.data[[cn]] <- k
}

# Clustree visualization across resolutions
ct <- clustree(seu@meta.data, prefix = "integrated_snn_res.", return = "plot")
ggsave(file.path(FIG_DIR, "Step1_clustree.png"), ct, width = 180, height = 140, units = "mm", dpi = 150)

# Choose a working resolution for downstream visualization/exports
RES_CHOSEN <- as.numeric(Sys.getenv("RES_CHOSEN", unset = "0.7"))
Idents(seu) <- seu$seurat_clusters <- factor(seu@meta.data[[paste0("integrated_snn_res.", RES_CHOSEN)]])

# UMAP split by sample, colored by chosen clusters
umap_split <- try({
  DimPlot(seu, reduction = "umap", split.by = "group",
          group.by = paste0("integrated_snn_res.", RES_CHOSEN), label = TRUE) +
    coord_fixed(ratio = 1)
}, silent = TRUE)
if (!inherits(umap_split, "try-error")) {
  ggsave(file.path(FIG_DIR, sprintf("Step1_umap_split_res%s.png", RES_CHOSEN)), umap_split,
         width = 180, height = 140, units = "mm", dpi = 150, bg = "white")
}

# QC violin plots (if the relevant columns exist)
qc_feats <- c("sum", "detected", "subsets_Mito_percent", "subsets_Rp_percent", "subsets_Heatshock_percent")
qc_feats <- qc_feats[qc_feats %in% colnames(seu@meta.data)]
if (length(qc_feats)) {
  vln <- VlnPlot(seu, features = qc_feats, pt.size = 0)
  ggsave(file.path(FIG_DIR, "Step1_qc_vln.png"), vln, width = 200, height = 140, units = "mm", dpi = 150, bg = "white")
}

discardCl <- seu$seurat_clusters %in% c(11, 12)
seu_filt <- seu[, !(discardCl)]
seu <- seu_filt

seu <- RunPCA(seu, npcs = 50, verbose = TRUE, features = PCA_features)

# Elbow (save to file for reproducibility)
elbow_plot <- ElbowPlot(seu, ndims = 50)
ggsave(file.path(FIG_DIR, "Step1_elbow_plot.png"), elbow_plot, width = 140, height = 120, units = "mm", dpi = 150)

# UMAP/Neighbors (dims 1:20 as in your script)
plan(PLAN_STRATEGY, workers = NWORKERS)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:20)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)

seu <- FindClusters(seu, resolution = c(seq(0.1, 1.0, 0.1)))
# Relabel clusters from 1..K (avoid "0") per original practice
meta_cols <- grep("integrated_snn_res.", colnames(seu@meta.data), value = TRUE)
for (cn in meta_cols) {
  k <- seu@meta.data[[cn]]
  levels(k) <- as.character(seq_len(nlevels(k)))
  seu@meta.data[[cn]] <- k
}

# Clustree visualization across resolutions
ct <- clustree(seu@meta.data, prefix = "integrated_snn_res.", return = "plot")
ggsave(file.path(FIG_DIR, "Step1_clustree.png"), ct, width = 180, height = 140, units = "mm", dpi = 150)

# Choose a working resolution for downstream visualization/exports
RES_CHOSEN <- as.numeric(Sys.getenv("RES_CHOSEN", unset = "0.7"))
Idents(seu) <- seu$seurat_clusters <- factor(seu@meta.data[[paste0("integrated_snn_res.", RES_CHOSEN)]])

# UMAP split by sample, colored by chosen clusters
umap_split <- try({
  DimPlot(seu, reduction = "umap", split.by = "group",
          group.by = paste0("integrated_snn_res.", RES_CHOSEN), label = TRUE) +
    coord_fixed(ratio = 1)
}, silent = TRUE)
if (!inherits(umap_split, "try-error")) {
  ggsave(file.path(FIG_DIR, sprintf("Step1_umap_split_res%s.png", RES_CHOSEN)), umap_split,
         width = 180, height = 140, units = "mm", dpi = 150, bg = "white")
}

# QC violin plots (if the relevant columns exist)
qc_feats <- c("sum", "detected", "subsets_Mito_percent", "subsets_Rp_percent", "subsets_Heatshock_percent")
qc_feats <- qc_feats[qc_feats %in% colnames(seu@meta.data)]
if (length(qc_feats)) {
  vln <- VlnPlot(seu, features = qc_feats, pt.size = 0)
  ggsave(file.path(FIG_DIR, "Step1_qc_vln.png"), vln, width = 200, height = 140, units = "mm", dpi = 150, bg = "white")
}

# Save integrated object (date-less by default)
res_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("%s_%s.rds", RESULTS_BASENAME, format(Sys.Date(), "%Y%m%d"))
            else sprintf("%s.rds", RESULTS_BASENAME)
saveRDS(seu, file = file.path(RES_DIR, res_file))

# -------------------- Tabulate cluster composition ------------------------
# Cell type label column may vary (e.g., ImmGen or Azimuth). Try common fields.
ctype_col <- intersect(c("CellType_immgen", "predicted.celltype", "predicted.annotation"), colnames(seu@meta.data))[1]

cluster_annot <- NULL
if (!is.na(ctype_col)) {
  cluster_annot <- tibble(Cluster = seu$seurat_clusters, Source = seu$group, CellType = seu@meta.data[[ctype_col]]) %>%
    group_by(Cluster, CellType) %>% summarise(no.cell = n(), .groups = "drop_last") %>%
    group_by(Cluster) %>% mutate(total.no = sum(no.cell), perc = 100 * no.cell / total.no) %>%
    arrange(Cluster, dplyr::desc(perc)) %>%
    top_n(n = 5, wt = perc)
}

source_cluster <- tibble(Cluster = seu$seurat_clusters, Group = seu$group) %>%
  group_by(Cluster, Group) %>% summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Group) %>% mutate(total.no = sum(no.cell), perc = 100 * no.cell / total.no) %>%
  ungroup() %>% select(Cluster, Group, perc)

# CC phase distribution (if present)
source_CCphase <- NULL
if ("CCphase" %in% colnames(seu@meta.data)) {
  source_CCphase <- tibble(Cluster = seu$seurat_clusters, CCphase = seu$CCphase) %>%
    group_by(Cluster, CCphase) %>% summarise(no.cell = n(), .groups = "drop_last") %>%
    group_by(CCphase) %>% mutate(total.no = sum(no.cell), perc = 100 * no.cell / total.no) %>%
    ungroup() %>% select(Cluster, CCphase, perc)
}

# Composition barplot (optional)
comp_plot <- ggplot(source_cluster, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") + coord_fixed(ratio = 1/10) + theme_bw() + xlab("group") + ylab("%")
ggsave(file.path(FIG_DIR, "Step1_cluster_composition.png"), comp_plot, width = 160, height = 120, units = "mm", dpi = 150)

# -------------------- Export summary tables (xlsx) ------------------------
wb <- wb_workbook()
if (!is.null(cluster_annot)) {
  wb <- wb_add_worksheet(wb, sheet = "cluster_annot")
  wb <- wb_add_data(wb, sheet = "cluster_annot", x = as.data.frame(cluster_annot))
}
source_cluster_wide <- source_cluster %>% tidyr::pivot_wider(names_from = Cluster, values_from = perc)
wb <- wb_add_worksheet(wb, sheet = "source_cluster")
wb <- wb_add_data(wb, sheet = "source_cluster", x = as.data.frame(source_cluster_wide))
if (!is.null(source_CCphase)) {
  source_CCphase_wide <- source_CCphase %>% tidyr::pivot_wider(names_from = Cluster, values_from = perc)
  wb <- wb_add_worksheet(wb, sheet = "source_CCphase")
  wb <- wb_add_data(wb, sheet = "source_CCphase", x = as.data.frame(source_CCphase_wide))
}
wb_save(wb, file.path(RES_DIR, "Step1_cluster_summaries.xlsx"), overwrite = TRUE)

# -------------------- Save final object & session info ---------------------
saveRDS(seu, file = file.path(RES_DIR, res_file))
writeLines(capture.output(sessionInfo()), file.path(RES_DIR, if (INCLUDE_DATE_IN_NAMES) sprintf("session_info_%s.txt", format(Sys.Date(), "%Y%m%d")) else "session_info.txt"))
