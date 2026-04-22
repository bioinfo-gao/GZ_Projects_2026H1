# Bioinformatics Analysis Report

Date: 2026-04-22 Project: RNAseq

## 1. Overview

This report summarizes the differential expression analysis and quality control metrics for the RNA-seq dataset. - **Analysis Tool**: DESeq2 - **Normalization**: VST (Variance Stabilizing Transformation) for PCA/Heatmap, Median-of-ratios for DE - **Significance Thresholds**: padj \< 0.05, \|log2FoldChange\| \>= 1

## 2. Quality Control (QC)

-   QC reports were generated using MultiQC.

## 3. Differential Expression Analysis Results

### Contrast: MQ-07-99_vs_DMSO

-   Total Significant Genes: 2970
    -   Upregulated: 1408
    -   Downregulated: 1562
-   Output File: `DEG_MQ-07-99_vs_DMSO.csv`

### Contrast: VK-8-101_vs_DMSO

-   Total Significant Genes: 2012
    -   Upregulated: 1346
    -   Downregulated: 666
-   Output File: `DEG_VK-8-101_vs_DMSO.csv`

### Contrast: MQ-07-81_vs_DMSO

-   Total Significant Genes: 1
    -   Upregulated: 1
    -   Downregulated: 0
-   Output File: `DEG_MQ-07-81_vs_DMSO.csv`

### Contrast: VK-8-101_vs_MQ-07-99

-   Total Significant Genes: 3033
    -   Upregulated: 1756
    -   Downregulated: 1277
-   Output File: `DEG_VK-8-101_vs_MQ-07-99.csv`

### Contrast: MQ-07-81_vs_MQ-07-99

-   Total Significant Genes: 2905
    -   Upregulated: 1537
    -   Downregulated: 1368
-   Output File: `DEG_MQ-07-81_vs_MQ-07-99.csv`

### Contrast: MQ-07-81_vs_VK-8-101

-   Total Significant Genes: 2066
    -   Upregulated: 683
    -   Downregulated: 1383
-   Output File: `DEG_MQ-07-81_vs_VK-8-101.csv`

## 4. Visualizations

### Principal Component Analysis (PCA)

-   **File**: `PCA.pdf`
-   **Description**: Shows sample clustering based on the top 500 most variable genes. Samples should cluster by biological group if the treatment effect is strong.

### Volcano Plots

-   **Files**: `Volcano_*.png`
-   **Description**: Displays the relationship between statistical significance (-log10 padj) and magnitude of change (log2FC). Red points indicate upregulated genes, blue points indicate downregulated genes.

## 5. Generated Data Files

| File Name | Description |
|:-----------------------------------|:-----------------------------------|
| `All_sample_gene_counts.tsv` | Raw count matrix for all samples. |
| `All_sample_gene_tpm.tsv` | TPM (Transcripts Per Million) matrix, if available. |
| `Normalized_Counts.csv` | DESeq2 normalized counts for downstream analysis. |
| `DEG_*.csv` | Detailed differential expression results including log2FC, p-values, and base means. |
| `PCA.pdf` | PCA plot showing sample relationships. |
| `Volcano_*.png` | Volcano plots for each contrast. |
| `QC/` | Directory containing MultiQC and other QC reports. |