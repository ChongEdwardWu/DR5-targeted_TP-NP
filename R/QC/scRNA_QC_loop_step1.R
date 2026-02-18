#!/usr/bin/env Rscript
# Minimal, looped QC pipeline for 10x Genomics scRNA-seq (Cell Ranger outputs).

suppressPackageStartupMessages({
  library(DropletUtils)
  library(BiocParallel)
  library(scran)
  library(scater)
  library(ggplot2)
  library(gridExtra)
  library(scDblFinder)
  library(bluster)
  library(BiocSingular)
})

set.seed(123)

## ---------------------- User-configurable parameters ----------------------
# You can override these via environment variables, e.g.:
#   export PROJECT_BASE="/path/to/project"
#   export NWORKERS=16
#   export INCLUDE_DATE_IN_LOG=FALSE
#   export NOTIFY_CMD="mail -s 'QC done' you@example.com <<< 'Finished.'"

PROJECT_BASE         <- Sys.getenv("PROJECT_BASE", unset = ".")        # Project root (default: current dir)
CRDIR_NAME           <- Sys.getenv("CRDIR_NAME", unset = "02_CRcount") # Subdir containing Cell Ranger outputs
OUTDIR_NAME          <- Sys.getenv("OUTDIR_NAME", unset = "03_R/QC")   # Where to write per-sample results
SAMPLE_FILTER_REGEX  <- Sys.getenv("SAMPLE_FILTER_REGEX", unset = ".*")# Optional regex to subset samples
NWORKERS             <- as.integer(Sys.getenv("NWORKERS", unset = "8"))
INCLUDE_DATE_IN_LOG  <- as.logical(Sys.getenv("INCLUDE_DATE_IN_LOG", unset = "FALSE"))
NOTIFY_CMD           <- Sys.getenv("NOTIFY_CMD", unset = "")           # Optional shell command for notification

# Optional: if you keep a personal startup file, point to a generic path here (left commented intentionally)
# source("~/.radian_profile")

# Parallel backend
bp <- MulticoreParam(workers = NWORKERS)

## ---------------------- Discover samples & loop ---------------------------
# Expecting Cell Ranger folders like: <PROJECT_BASE>/<CRDIR_NAME>/<sample>/outs/filtered_feature_bc_matrix
Basedir <- normalizePath(PROJECT_BASE, winslash = "/", mustWork = FALSE)
CRdir   <- file.path(Basedir, CRDIR_NAME)

# List all immediate subdirectories in CRdir as samples, with optional regex filter
samples <- list.files(CRdir, pattern = SAMPLE_FILTER_REGEX, full.names = FALSE)
samples <- samples[sapply(file.path(CRdir, samples), dir.exists)]

for (sample in samples) {
  message("========== Analyzing sample: ", sample, " ==========")

  ## ---------------------- Prepare per-sample folders & logging ------------
  workdir <- file.path(Basedir, OUTDIR_NAME, sample)
  figdir  <- file.path(workdir, "figures")
  if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)
  if (!dir.exists(figdir))  dir.create(figdir)

  # Log file: by default, do NOT include a date to avoid revealing run time in filenames
  log_stub <- sprintf("%s_scRNA_QC_step1", sample)
  log_file <- if (isTRUE(INCLUDE_DATE_IN_LOG)) {
    file.path(workdir, sprintf("%s_%s.log", log_stub, format(Sys.time(), "%Y%m%d")))
  } else {
    file.path(workdir, sprintf("%s.log", log_stub))
  }

  # Redirect console output/messages to the log; keeps repo clean and reproducible
  output_connection <- file(log_file, open = "wt")
  sink(output_connection)
  sink(output_connection, type = "message")

  cat(sprintf("[%s] Starting QC for %s\n", as.character(Sys.time()), sample))

  ## ---------------------- Input: 10x count matrix -------------------------
  # Standard Cell Ranger location for filtered counts:
  inpath <- file.path(CRdir, sample, "outs", "filtered_feature_bc_matrix")
  message("Input 10X path: ", inpath)

  # Read counts into a SingleCellExperiment
  sce <- read10xCounts(inpath, col.names = TRUE)

  ## ---------------------- Empty droplets (cell calling) --------------------
  # Use the 2nd-smallest total UMI as a conservative lower bound (common heuristic)
  a <- colSums(counts(sce))
  lower_bound <- sort(a)[2]
  e.out <- emptyDrops(counts(sce), lower = lower_bound, BPPARAM = bp)
  print(summary(e.out$FDR <= 0.05))

  # QC plot: emptyDrops significance vs total
  png(file = file.path(figdir, "01_Cell_Calling_for_Droplets.png"), width = 250, height = 250, units = "mm", res = 150)
  plot(e.out$Total, -e.out$LogProb,
       col = ifelse(e.out$FDR <= 0.05, "red", "black"),
       xlab = "Total UMI count", ylab = "-Log Probability")
  dev.off()

  print(table(Sig = e.out$FDR <= 0.05, Limited = e.out$Limited))

  # Keep only detected cells at FDR <= 0.05 (tunable)
  sce <- sce[, which(e.out$FDR <= 0.05)]

  ## ---------------------- Normalization & HVGs -----------------------------
  # Quick clustering for deconvolution-based size factors
  clust <- quickCluster(sce, BPPARAM = bp)
  print(table(clust))

  # (Optional inspection) deconvolution size factors; not used further below.
  # NOTE: We then call computeSumFactors() for normalization and log-transform.
  deconv.sf <- calculateSumFactors(sce, cluster = clust, BPPARAM = bp)
  print(summary(deconv.sf))

  # Library size normalization + log transform
  norm <- computeSumFactors(sce, cluster = clust, min.mean = 0.1, BPPARAM = bp)
  norm <- logNormCounts(norm)

  # Model gene variance, pick top HVGs (n=2000 is a common default)
  dec <- modelGeneVar(norm)
  fit <- metadata(dec)

  png(file = file.path(figdir, "02_Feature_selection.png"), width = 250, height = 250, units = "mm", res = 150)
  par(mfrow = c(1, 1))
  plot(fit$mean, fit$var, xlab = "Mean of log-expression", ylab = "Variance of log-expression")
  curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
  dev.off()

  hvg.var <- getTopHVGs(dec, n = 2000)
  print(length(hvg.var))

  ## ---------------------- Dimensionality reduction & clustering -----------
  # PCA on HVGs, using a randomized SVD for speed on large matrices
  norm <- fixedPCA(norm, subset.row = hvg.var, BSPARAM = RandomParam(), name = "PCA")

  # UMAP projection
  norm <- runUMAP(norm, dimred = "PCA", BPPARAM = bp)

  # Graph-based clustering on PCA space (Jaccard NN graph)
  colLabels(norm) <- clusterCells(norm, use.dimred = "PCA", BLUSPARAM = NNGraphParam(type = "jaccard"))

  # UMAP plot with cluster labels
  sceumap <- plotReducedDim(norm, "UMAP", colour_by = "label", text_by = "label")
  ggsave(filename = file.path(figdir, "03_Dimensionality_Reduction_and_Clustering.png"),
         plot = sceumap, width = 250, height = 250, units = "mm", dpi = 150, device = "png", bg = "white")

  ## ---------------------- Doublet detection --------------------------------
  # Use scDblFinder with pre-clustering info to stabilize estimates
  dbl <- scDblFinder(norm, clusters = clust, nfeatures = 2000, BPPARAM = bp)
  print(table(dbl$scDblFinder.class))
  norm$DoubletScore <- dbl$scDblFinder.score
  norm$DoubletClass <- dbl$scDblFinder.class

  ## ---------------------- Basic QC metrics ---------------------------------
  # Define gene sets (regex covers human/mouse mitochondrial, ribosomal, heat-shock, hemoglobin, platelet markers)
  is.mito       <- grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l", rowData(norm)$Symbol)
  is.rp         <- grep("^Rp[sl]|^RP[SL]", rowData(norm)$Symbol)
  is.hsp        <- grep("^HSP|^DNAJ|^Hsp|^Dnaj", rowData(norm)$Symbol)
  is.hemoglobin <- grep("^Hb[ab]-|^HB[AB]", rowData(norm)$Symbol)
  is.plat       <- grep("Pecam1|Pf4|PECAM1|PF4", rowData(norm)$Symbol)

  qc <- perCellQCMetrics(norm, subsets = list(
    Mito = is.mito, Heatshock = is.hsp, Rp = is.rp, RBC = is.hemoglobin, Plat = is.plat
  ))

  # Attach QC metrics to colData
  colData(norm) <- cbind(colData(norm), qc)

  # Fixed-threshold filters (adjust to your dataset/species as needed)
  qc.nexprs <- qc$detected < 100
  qc.mito   <- qc$subsets_Mito_percent > 15
  qc.rbc    <- qc$subsets_RBC_percent  > 1
  qc.plat   <- qc$subsets_Plat_percent > 1
  fix.discard <- qc.nexprs | qc.mito | qc.rbc | qc.plat
  print(table(fix.discard))
  norm$fix.discard <- fix.discard

  # Mito vs detected features (flagged cells highlighted)
  png(file = file.path(figdir, "04_Cell_Filtering.png"), width = 250, height = 250, units = "mm", res = 150)
  plot(qc$detected, qc$subsets_Mito_percent, log = "x",
       xlab = "Detected features", ylab = "Mitochondrial %")
  points(qc$detected[fix.discard], qc$subsets_Mito_percent[fix.discard], col = "dodgerblue", pch = 16, cex = 0.5)
  legend("topright", legend = c("Fix", "Adaptive"), col = c("dodgerblue", "red"), lty = 2, cex = 1)
  dev.off()

  # Final discard flag (here equal to fixed thresholds; extend with adaptive rules if desired)
  discard <- fix.discard
  norm$discard <- discard

  # Diagnostic panels across QC metrics
  diagnostic <- gridExtra::grid.arrange(
    plotColData(norm, x = "DoubletClass", y = "sum", colour_by = "DoubletClass") + ggtitle("Doublets"),
    plotColData(norm, x = "discard", y = "sum",  colour_by = "discard") + ggtitle("Total count"),
    plotColData(norm, x = "discard", y = "subsets_Mito_percent", colour_by = "discard") + ggtitle("Mito percent"),
    plotColData(norm, x = "discard", y = "subsets_Rp_percent",   colour_by = "discard") + ggtitle("Ribosomal protein percent"),
    plotColData(norm, x = "discard", y = "subsets_Heatshock_percent", colour_by = "discard") + ggtitle("Heatshock protein percent"),
    plotColData(norm, x = "sum", y = "detected", colour_by = "discard") + ggtitle("Detected features"),
    plotColData(norm, x = "sum", y = "subsets_RBC_percent", colour_by = "discard") + ggtitle("Hemoglobin percent"),
    plotColData(norm, x = "sum", y = "subsets_Plat_detected", colour_by = "discard") + ggtitle("Platelet percent"),
    ncol = 3
  )
  ggsave(filename = file.path(figdir, "05_Diagnostic_Plots.png"),
         plot = diagnostic, width = 500, height = 250, units = "mm", dpi = 150, device = "png", bg = "white")

  # Per-cluster QC summaries
  clusQC <- gridExtra::grid.arrange(
    plotColData(norm, x = "label", y = "DoubletScore", colour_by = "DoubletClass") + ggtitle("Doublets"),
    plotColData(norm, x = "label", y = "sum",          colour_by = "label") + ggtitle("Total count"),
    plotColData(norm, x = "label", y = "detected",     colour_by = "label") + ggtitle("Detected features"),
    plotColData(norm, x = "label", y = "subsets_Mito_percent", colour_by = "label") + ggtitle("Mito percent"),
    plotColData(norm, x = "label", y = "subsets_Rp_percent",   colour_by = "label") + ggtitle("Ribosomal genes percent"),
    plotColData(norm, x = "label", y = "subsets_Heatshock_percent", colour_by = "label") + ggtitle("Heatshock protein genes percent"),
    plotColData(norm, x = "label", y = "subsets_RBC_percent",  colour_by = "label") + ggtitle("RBC percent"),
    plotColData(norm, x = "label", y = "subsets_Plat_percent", colour_by = "label") + ggtitle("Platelets percent"),
    plotColData(norm, x = "label", y = "discard",      colour_by = "label") + ggtitle("Discard"),
    ncol = 3
  )
  ggsave(filename = file.path(figdir, "06_Cluster_QC_Plots.png"),
         plot = clusQC, width = 500, height = 500, units = "mm", dpi = 150, device = "png", bg = "white")

  # UMAP views: kept vs discarded
  sceumap2 <- gridExtra::grid.arrange(
    plotReducedDim(norm[, norm$discard == 0], "UMAP", colour_by = "label", text_by = "label") + ggtitle("Cells remained"),
    plotReducedDim(norm[, norm$discard == 1], "UMAP", colour_by = "label", text_by = "label") + ggtitle("Cells discarded"),
    ncol = 2
  )
  ggsave(filename = file.path(figdir, "07_PreFiltering_UMAP.png"),
         plot = sceumap2, width = 250, height = 250, units = "mm", dpi = 150, device = "png", bg = "white")

  ## ---------------------- Save results -------------------------------------
  saveRDS(norm, file = file.path(workdir, sprintf("01_QC_step1_%s.rds", sample)))
  message(sprintf("Job %s completed at %s.", sample, as.character(Sys.time())))

  ## ---------------------- Cleanup & next sample ----------------------------
  sink(); sink(type = "message")
  rm(list = setdiff(ls(), c("Basedir", "CRdir", "NWORKERS", "bp", "samples",
                            "OUTDIR_NAME", "SAMPLE_FILTER_REGEX",
                            "INCLUDE_DATE_IN_LOG", "NOTIFY_CMD")))
  gc()
  graphics.off()
}

# End of Script
