# Bioinformatics Analysis Report

Date: 2026-05-03
Project: 2026_Item11_Mingxue

## 1. Overview
This report summarizes the differential expression analysis and quality control metrics for the RNA-seq dataset.
- **Analysis Mode**: Exploratory analysis without biological replicates
- **Analysis Tool**: DESeq2 size-factor normalization only; no formal hypothesis testing
- **Normalization**: Median-of-ratios normalization, followed by log2(normalized count + 1) for PCA/Heatmap
- **Highlight Thresholds**: |log2FoldChange| >= 1 (exploratory only; p-values/padj are not reported)

## 2. Quality Control (QC)
- QC reports were generated using MultiQC.

## 3. Differential Expression Analysis Results

### Contrast: Minxue_T_vs_Minxue_C
- Total Highlighted Genes: 630
  - Upregulated: 185
  - Downregulated: 445
- Output File: `DEG_Minxue_T_vs_Minxue_C.csv`

## 4. Visualizations

### Principal Component Analysis (PCA)
- **File**: `PCA.pdf`
- **Description**: Shows sample clustering based on the top variable genes used for PCA.

### Fold-change Plot
- **Files**: `Exploratory_*.png`
- **Description**: Exploratory scatter plot showing log2FC against mean normalized expression.

### Heatmap
- **File**: `heatmap.pdf`
- **Description**: Hierarchical clustering of the top 50 highest-change genes across all samples.

## 5. Generated Data Files

| File Name | Description |
| :--- | :--- |
| `All_sample_gene_counts.tsv` | Raw count matrix for all samples. |
| `All_sample_gene_tpm.tsv` | TPM matrix, if available. |
| `Normalized_Counts.csv` | DESeq2 normalized counts for downstream analysis. |
| `DEG_*.csv` | Exploratory fold-change table including normalized means, log2FC, and base means. |
| `PCA.pdf` | PCA plot showing sample relationships. |
| `Exploratory_*.png` | Exploratory fold-change scatter plots for each contrast. |
| `heatmap.pdf` | Heatmap of top genes. |
| `QC/` | Directory containing MultiQC and other QC reports. |
