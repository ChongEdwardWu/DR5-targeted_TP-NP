#!/usr/bin/env Rscript
# scRNA-seq — Step 3: Figures & summary exports

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(patchwork)
  library(tidyverse)
  library(ggplot2)
  library(clustree)
  library(openxlsx2)
  library(future)
  library(viridis)
  library(Nebulosa)
  library(UCell)
  library(ggpubr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(fgsea)
  library(msigdbr)
})

set.seed(123)
# ------------------------- Config -----------------------------------------
NWORKERS              <- as.integer(Sys.getenv("NWORKERS", unset = "8"))
PLAN_STRATEGY         <- Sys.getenv("PLAN", unset = "sequential")  # or "multicore"
PROJECT_BASE          <- Sys.getenv("PROJECT_BASE", unset = ".")
WORKDIR               <- Sys.getenv("WORKDIR", unset = file.path(PROJECT_BASE, "03_R"))
RES_DIR               <- file.path(WORKDIR, "results")
FIG_DIR               <- file.path(WORKDIR, "figures")
ANNOTATED_RDS         <- Sys.getenv("ANNOTATED_RDS", unset = file.path(RES_DIR, "02_Annotation.rds"))
INCLUDE_DATE_IN_NAMES <- as.logical(Sys.getenv("INCLUDE_DATE_IN_NAMES", unset = "FALSE"))
SPECIES               <- Sys.getenv("SPECIES", unset = "mm")
GROUP_LEVELS_ENV      <- Sys.getenv("GROUP_LEVELS", unset = "PBS,MD5,TP")

if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
setwd(WORKDIR)
plan(PLAN_STRATEGY, workers = NWORKERS)

# ------------------------- Load object ------------------------------------
seu <- readRDS(ANNOTATED_RDS)

# Consistent grouping (if available)
if ("group" %in% colnames(seu@meta.data)) {
  group_levels <- trimws(strsplit(GROUP_LEVELS_ENV, ",")[[1]])
  seu$group <- factor(seu$group, levels = group_levels)
}

# Workbook to collect summary tables
wb <- wb_workbook()

# ------------------------- Fig 1: UMAP --------------------------------
p1 <- DimPlot(seu, reduction = "umap", group.by = "CellType", label = FALSE) +
  coord_fixed(ratio = 1) + scale_color_manual(values = colorpanel)

if (!inherits(p1, "try-error")) {
  ggsave(file.path(FIG_DIR,  "fig1_UMAP.png"), p1, width = 250, height = 150, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR,  "fig1_UMAP.pdf"), p1, width = 250, height = 150, units = "mm", device = "pdf", bg = "transparent")
}

# ------------------------- Fig 2: marker VlnPlot -----------------------
DefaultAssay(seu) <- "SCT"
feature_genes <- list(
  Tcell = c("Cd3e", "Cd8a", "Cd4", "Foxp3"),
  Bcell = c("Cd79a", "Cd19"),
  NK = c("Ncr1", "Klrb1c"),
  DC = c("Batf3", "Itgax")
)

gene_panel <- unique(unlist(feature_genes))
gene_panel <- intersect(gene_panel, rownames(seu))


p2 <- VlnPlot(seu, group.by = "CellType", assay = "SCT", features = gene_panel, stack = T, pt.size=0.1)

ggsave(file.path(FIG_DIR, "fig2_CellTypeFeat_vlnplot.png"), p2, width = 100, height = 300, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR, "fig2_CellTypeFeat_vlnplot.pdf"), p2, width = 100, height = 300, units = "mm", device = "pdf", bg = "transparent")

# ------------------------- Fig 3: cluster abundance -----------------------
source_cluster <- tibble(Cluster = seu$CellType, Group = seu$group) %>%
  group_by(Cluster, Group) %>% summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Group) %>% mutate(total.no = sum(no.cell), perc = 100 * no.cell / total.no)

wb <- wb_add_worksheet(wb, sheet = "group_cluster")
wb <- wb_add_data(wb, sheet = "group_cluster", x = as.data.frame(source_cluster))

p3_bar <- ggplot(source_cluster, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") + scale_fill_manual(values = colorpanel) +
  coord_fixed(ratio = 1/10) + theme_bw() + xlab("group") + ylab("%")

# Save abundance plots
for (nm in c("fig3_cluster_abundance.png", "fig3_cluster_abundance.pdf")) {
  ggsave(file.path(FIG_DIR, nm), p3_bar, width = 100, height = 100, units = "mm", dpi = 300, device = tools::file_ext(nm))
}

# ------------------------- Fig 4: signature scores ( DC & CD8T) ------------------------
pathwaysGO <- msigdbr("mouse", category = "C5")
pathwaysGO <- split(pathwaysGO$gene_symbol, pathwaysGO$gs_name)

pathwaysImm <- msigdbr("mouse", category = "C7")
pathwaysImm <- split(pathwaysImm$gene_symbol, pathwaysImm$gs_name)

# violin plot
grp_fill <- c(PBS = "#5F5F5F", MD5 = "#84c4b7", TP = "#ff8080")
make_sig_violin <- function(sig){
  dat <- dplyr::filter(df_long, Signature == sig)
  y_top <- max(dat$Score, na.rm = TRUE)
  y_rng <- diff(range(log(dat$Score)+1, na.rm = TRUE))
  y_lab <- y_top + 0.06 * ifelse(is.finite(y_rng) && y_rng > 0, y_rng, abs(y_top) + 1)
  ggplot(dat, aes(x = group, y = Score, fill = group)) +
    geom_violin(trim = FALSE, scale = "area", color = "black") +
    geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
    stat_summary(fun = median, geom = "point", size = 1) +
    scale_fill_manual(values = grp_fill) +
    labs(x = NULL, y = "UCell score", title = sig) +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"), aspect.ratio = 1.5) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) 
}

DefaultAssay(seu) <- "SCT"

if ("CellType" %in% colnames(seu@meta.data)) {
  cells.use <- WhichCells(seu, expression = CellType == "DC")
  seu_sub   <- subset(seu, cells = cells.use)
  if (ncol(seu_sub) > 0) {
    seu_sub$group <- factor(seu_sub$group)

# Define gene signatures
sig_list <- list(
  Phagocytosis_recognition   = pathwaysGO["GOBP_PHAGOCYTOSIS_RECOGNITION"],
  MHC_complex  = c("H2-Kd", "H2-Dd", "H2-Ld", "H2-Aa", "H2-Ab1", "H2-Ea", "H2-Eb1",
                                  "H2-DMa", "H2-DMb1", "Tap1", "Tap2", "Psmb8", "Psmb9", "B2m")
)

# Score (UCell). Columns appended as *_UCell
seu_sub <- AddModuleScore_UCell(seu_sub, features = sig_list, ncores = NWORKERS, maxRank = 1500)

# vlnplot per signature
sig_for_violin <- grep("_UCell$", colnames(seu_sub@meta.data), value = TRUE)

df_long <- FetchData(seu_sub, vars = c("group", sig_for_violin)) %>%
  as_tibble() %>%
  pivot_longer(all_of(sig_for_violin), names_to = "Signature", values_to = "Score") %>%
  mutate(group = factor(group, levels = c("PBS", "MD5", "TP")),
         Signature = factor(Signature, levels = sig_for_violin))

p_list <- lapply(sig_for_violin, make_sig_violin)
p_DC <- ggpubr::ggarrange(plotlist = p_list, ncol = 2, nrow = 1, align = "hv")

  ggsave(file.path(FIG_DIR, sprintf("fig4_DC_signature_vlnplot_%s.png", sig)), p_DC, width = 100, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR, sprintf("fig4_DC_signature_vlnplot_%s.pdf", sig)), p_DC, width = 100, height = 100, units = "mm", device = "pdf", bg = 'transparent')
 }
}

if ("CellType" %in% colnames(seu@meta.data)) {
  cells.use <- WhichCells(seu, expression = CellType == "CD8T")
  seu_sub   <- subset(seu, cells = cells.use)
  if (ncol(seu_sub) > 0) {
    seu_sub$group <- factor(seu_sub$group)

# Define gene signatures
sig_list <- list(
  Tcm   = c("Igfbp4","Sell", "Dapl1", "Ccr7","Lef1", "Tcf7", "S1pr1", "Ccl5", "Id2"),
  TeffvsTex  = pathwaysImm["GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP"]
)

# Score (UCell). Columns appended as *_UCell
seu_sub <- AddModuleScore_UCell(seu_sub, features = sig_list, ncores = NWORKERS, maxRank = 1500)

# vlnplot per signature
sig_for_violin <- grep("_UCell$", colnames(seu_sub@meta.data), value = TRUE)

df_long <- FetchData(seu_sub, vars = c("group", sig_for_violin)) %>%
  as_tibble() %>%
  pivot_longer(all_of(sig_for_violin), names_to = "Signature", values_to = "Score") %>%
  mutate(group = factor(group, levels = c("PBS", "MD5", "TP")),
         Signature = factor(Signature, levels = sig_for_violin))

p_list <- lapply(sig_for_violin, make_sig_violin)
p_CD8T <- ggpubr::ggarrange(plotlist = p_list, ncol = 2, nrow = 1, align = "hv")

  ggsave(file.path(FIG_DIR, sprintf("fig4_CD8T_signature_vlnplot_%s.png", sig)), p_CD8T, width = 100, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR, sprintf("fig4_CD8T_signature_vlnplot_%s.pdf", sig)), p_CD8T, width = 100, height = 100, units = "mm", device = "pdf", bg = 'transparent')
 }
}

# ------------------------- Fig 5: TP-NPs vs MD5-1 GSEA ( DC & CD8T) ------------
    # Exclude confounders from DE
    rn <- if ("RNA" %in% names(seu_sub@assays)) rownames(seu_sub@assays$RNA@counts) else rownames(seu_sub)
    hist_genes    <- grep("^Hist", rn, ignore.case = TRUE, value = TRUE)
    hb_genes      <- grep("^Hb[ab]-|^HB[^(P)]", rn, value = TRUE)
    mt_genes      <- grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l", rn, ignore.case = TRUE, value = TRUE)
    rps_genes     <- grep("^Rp[sl]|^RP[SL]", rn, ignore.case = TRUE, value = TRUE)
    rik_genes     <- grep("^Rik", rn, ignore.case = TRUE, value = TRUE)
    alu_genes     <- grep("^AL", rn, ignore.case = TRUE, value = TRUE)
    pseudo_genes  <- grep("-rs|-ps", rn, ignore.case = TRUE, value = TRUE)
    mir_genes     <- grep("^Mir", rn, ignore.case = TRUE, value = TRUE)
    gencode_genes <- grep("^Gm", rn, ignore.case = TRUE, value = TRUE)

    bad_features <- unique(c(hist_genes, hb_genes, mt_genes, rps_genes, rik_genes, alu_genes, pseudo_genes, mir_genes, gencode_genes))
    features_de  <- setdiff(rn, bad_features)

if ("CellType" %in% colnames(seu@meta.data)) {
  cells.use <- WhichCells(seu, expression = CellType == "DC")
  seu_sub   <- subset(seu, cells = cells.use)
  if (ncol(seu_sub) > 0) {
    seu_sub$group <- factor(seu_sub$group)

    DefaultAssay(seu_sub) <- "RNA"
    seu_sub <- NormalizeData(seu_sub)
    Idents(seu_sub) <- "group"

    deg <- FindMarkers(seu_sub, assay = "RNA", features = features_de, ident.1 = "TP", ident.2 = "MD5",
                       test.use = "wilcox", min.pct = 0.10, logfc.threshold = 0.05, densify = TRUE) %>% 
                  as.data.frame() %>% rownames_to_column("Gene") %>%
		  arrange(desc(avg_log2FC))

   logFC <- deg$avg_log2FC
   names(logFC) <- deg[, 1]
   gsea_list <- sort(logFC, decreasing = TRUE)
   gsea_list <- gsea_list[!is.na(names(gsea_list))]

   gsea_result <- GSEA(geneList = gsea_list, exponent = 1, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.5, pAdjustMethod = "BH", TERM2GENE = pathwaysGO, verbose = TRUE)

geneset <- c( 
  "GOBP_INNATE_IMMUNE_RESPONSE",
  "GOCC_RECYCLING_ENDOSOME_MEMBRANE", 
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN", 
  "GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION"
)
p_DC <- gseaplot2( gsea_result,  geneSetID = geneset,  base_size = 12)

  ggsave(file.path(FIG_DIR, sprintf("fig5_DC_signature_GSEA_%s.png", sig)), p_DC, width = 100, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR, sprintf("fig5_DC_signature_GSEA_%s.pdf", sig)), p_DC, width = 100, height = 100, units = "mm", device = "pdf", bg = 'transparent')
  }
}

if ("CellType" %in% colnames(seu@meta.data)) {
  cells.use <- WhichCells(seu, expression = CellType == "CD8T")
  seu_sub   <- subset(seu, cells = cells.use)
  if (ncol(seu_sub) > 0) {
    seu_sub$group <- factor(seu_sub$group)

    DefaultAssay(seu_sub) <- "RNA"
    seu_sub <- NormalizeData(seu_sub)
    Idents(seu_sub) <- "group"

    deg <- FindMarkers(seu_sub, assay = "RNA", features = features_de, ident.1 = "TP", ident.2 = "MD5",
                       test.use = "wilcox", min.pct = 0.10, logfc.threshold = 0.05, densify = TRUE) %>% 
                  as.data.frame() %>% rownames_to_column("Gene") %>%
		  arrange(desc(avg_log2FC))

   logFC <- deg$avg_log2FC
   names(logFC) <- deg[, 1]
   gsea_list <- sort(logFC, decreasing = TRUE)
   gsea_list <- gsea_list[!is.na(names(gsea_list))]

   gsea_result <- GSEA(geneList = gsea_list, exponent = 1, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.5, pAdjustMethod = "BH", TERM2GENE = pathwaysGO, verbose = TRUE)

geneset <- c( 
  "GOBP_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION", 
  "GOBP_LEUKOCYTE_DEGRANULATION", 
  "GOBP_T_CELL_PROLIFERATION",
  "GOBP_LYMPHOCYTE_MIGRATION",
  "GOBP_IMMUNOLOGICAL_MEMORY_PROCESS",
  "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"
)
p_CD8T <- gseaplot2( gsea_result,  geneSetID = geneset,  base_size = 12)

  ggsave(file.path(FIG_DIR, sprintf("fig5_CD8T_signature_GSEA_%s.png", sig)), p_CD8T, width = 100, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR, sprintf("fig5_CD8T_signature_GSEA_%s.pdf", sig)), p_CD8T, width = 100, height = 100, units = "mm", device = "pdf", bg = 'transparent')
  }
}

# ------------------------- Save outputs -----------------------------------
wb_save(wb, file.path(RES_DIR, "Step3_results.xlsx"), overwrite = TRUE)

res_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("03_Figures%s.rds", format(Sys.Date(), "%Y%m%d")) else "03_Figures.rds"
saveRDS(seu, file = file.path(RES_DIR, res_file))

sess_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("session_info_step3_%s.txt", format(Sys.Date(), "%Y%m%d")) else "session_info_step3.txt"
writeLines(capture.output(sessionInfo()), file.path(RES_DIR, sess_file))
