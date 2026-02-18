# DR5-targeted_TP-NP

Code for the single-cell RNA-seq analysis used in the manuscript  
**“Single-component receptor-density-gated fibrillar scaffold assembly drives higher-order DR5 clustering to induce tumour-selective pyroptosis.”**

This repository contains R scripts that reproduce the core scRNA-seq processing and analysis workflow, including quality control, doublet removal, integration across treatment groups, clustering and annotation, differential expression and signature scoring, enrichment analyses, and figure generation for the manuscript.

---

## Repository structure

- `R/QC/scRNA_QC_loop_step1.R`  
  Sample-level preprocessing and quality control starting from Cell Ranger outputs. This script performs:
  - loading raw and filtered matrices,
  - `emptyDrops`-based cell calling,
  - QC metric calculation,
  - preliminary filtering,
  - doublet detection with `scDblFinder`,
  - initial clustering for QC inspection, and
  - saving per-sample intermediate Seurat objects.  
  Inputs (per sample): Cell Ranger `raw_feature_bc_matrix` and `filtered_feature_bc_matrix` under `02_CRcount/<SAMPLE>/outs/` (default).  
  Outputs (per sample, default): `03_R/QC/<SAMPLE>/01_QC_step1_<SAMPLE>.rds`

- `R/QC/scRNA_QC_loop_step2.R`  
  Per-sample final QC and Seurat object preparation for downstream integration. This script performs:
  - loading the Step 1 QC object,
  - final filtering,
  - cell-cycle scoring,
  - optional reference-style annotation helpers (e.g., `SingleR`, `scMCA`, depending on your setup),
  - `SCTransform` normalisation, and
  - saving the final per-sample Seurat object used for integration.  
  Input: `03_R/QC/<SAMPLE>/01_QC_step1_<SAMPLE>.rds`  
  Output (default): `03_R/QC/<SAMPLE>/<SAMPLE>_seu.rds`  
  Note: this script is intended to be run once per sample with `SAMPLE` specified.

- `R/Analysis/Step1_Integration.R`  
  Integration across samples and treatment groups (Seurat anchor-based integration). This script:
  - loads per-sample Seurat objects created by QC Step 2,
  - assigns group metadata (default expected sample IDs are `PBS_1–3`, `MD5_1–3`, and `TP_1–3`),
  - performs integration, dimensional reduction, UMAP, and initial clustering,
  - supports optional removal of low-quality clusters, and
  - writes an integrated Seurat object for downstream analysis.  
  Inputs: `03_R/QC/<SAMPLE>/<SAMPLE>_seu.rds`  
  Main output (default): `03_R/results/01_Integration.rds`

- `R/Analysis/Step2_Clustering.R`  
  Major-lineage clustering and annotation starting from the integrated object. This script:
  - refines clustering,
  - identifies marker genes,
  - annotates major immune lineages, and
  - produces core diagnostic plots.  
  Input: `03_R/results/01_Integration.rds`  
  Output (default): `03_R/results/02_Clustering.rds`

- `R/Analysis/Step3_Figures.R`  
  Downstream analysis and figure-generation code used for manuscript figures, including:
  - module scoring and signature scoring (e.g., `UCell`),
  - differential expression in selected comparisons,
  - pathway enrichment and GSEA-style analyses, and
  - generation of publication plots.  
  Inputs: the annotated Seurat object from Step 2 and additional intermediate objects created within the script.  
  Note: this script contains several analysis- and plot-specific settings (for example, figure output paths, palettes, and gene panels) and may require minor edits to run cleanly in a new environment.

---

## Analysis workflow

1. **Prepare Cell Ranger outputs**  
   Place Cell Ranger outputs for each sample under the repository root using the default layout:
   - `02_CRcount/<SAMPLE>/outs/raw_feature_bc_matrix/`
   - `02_CRcount/<SAMPLE>/outs/filtered_feature_bc_matrix/`

2. **Run QC Step 1**  
   From the repository root:
   ```bash
   Rscript R/QC/scRNA_QC_loop_step1.R
   ```
   This generates per-sample QC objects under `03_R/QC/<SAMPLE>/`.

3. **Run QC Step 2 for each sample**  
   Example:
   ```bash
   SAMPLE=PBS_1 Rscript R/QC/scRNA_QC_loop_step2.R
   ```
   Repeat for all samples (by default the scripts expect `PBS_1–3`, `MD5_1–3`, and `TP_1–3`, but you can adapt the sample list).

4. **Integrate samples**  
   ```bash
   Rscript R/Analysis/Step1_Integration.R
   ```
   Output: `03_R/results/01_Integration.rds`

5. **Clustering and annotation**  
   ```bash
   Rscript R/Analysis/Step2_Clustering.R
   ```
   Output: `03_R/results/02_Clustering.rds`

6. **Downstream analyses and figure generation**  
   ```bash
   Rscript R/Analysis/Step3_Figures.R
   ```
   This script produces analysis outputs and figures used in the manuscript. It may require adjusting a small number of script-level parameters (paths and plot settings) depending on your environment.

---

## Input data

This repository does not include raw sequencing data or large intermediate objects.

To reproduce the workflow, you will need:
- Cell Ranger outputs (raw and filtered matrices) for each sample, organised as described above, or equivalent matrices adapted to the scripts.
- Sufficient compute and memory for `SCTransform`, integration, and downstream plotting.

Please refer to the manuscript for details on the experimental design, sample collection, and data availability. If accession numbers are provided in the final published version, they should be used as the canonical source for raw data retrieval.

---

## Configuration and paths

Several scripts use environment variables to simplify portability. Defaults are defined in each script.

Common variables include:
- `PROJECT_BASE`  
  Project root directory. Defaults to the current working directory.
- `CRDIR_NAME`  
  Directory name that contains Cell Ranger outputs. Default is `02_CRcount`.
- `WORKDIR`  
  Working directory for outputs. Default is `PROJECT_BASE/03_R`.
- `QC_DIR`, `RES_DIR`, `FIG_DIR`  
  Output subdirectories under `WORKDIR`.
- `SAMPLE`  
  Required for `scRNA_QC_loop_step2.R` to select the sample to process.

If your data are stored in a different layout, edit the path definitions near the top of each script.

---

## Software requirements

The analysis was developed in **R** (recent 4.x series). Key packages used across scripts include (non-exhaustive):
- `Seurat`, `SeuratObject`, `sctransform`
- `DropletUtils`, `scater`, `scran`
- `scDblFinder`
- `SingleR` (optional), `scMCA` (optional), `Azimuth` (optional, depending on your setup)
- `dplyr`, `data.table`, `Matrix`, `stringr`, `tidyr`
- `ggplot2`, `patchwork`, `cowplot`
- `UCell`
- `clusterProfiler`, `fgsea`, `msigdbr`, `org.Mm.eg.db` (for enrichment analyses)

Please consult the `library()` calls in each script for the exact dependencies.

---

## Citation

If you use this code, please cite the associated manuscript:

**Single-component receptor-density-gated fibrillar scaffold assembly drives higher-order DR5 clustering to induce tumour-selective pyroptosis.**  
manuscript submitted.

---

## Contact

For questions about the code or analysis, please open an issue in this GitHub repository or contact:

- **Chong Wu** – wuchong5@mail(dot)sysu(dot)edu(dot)cn
