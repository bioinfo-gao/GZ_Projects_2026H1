#!/usr/bin/env Rscript

# Load required libraries
library(DESeq2)
library(dplyr)
library(readr)

# Set working directory
setwd("/home/gao/projects/tmp")

# Read the existing DEG files to extract raw counts
df_g2_g1 <- read_csv("DEG_results_group2-vs-group1.csv")
df_g3_g1 <- read_csv("DEG_results_group3-vs-group1.csv") 
df_g4_g1 <- read_csv("DEG_results_group4-vs-group1.csv")

# Set Geneid as rownames for proper alignment
rownames(df_g2_g1) <- df_g2_g1$Geneid
rownames(df_g3_g1) <- df_g3_g1$Geneid  
rownames(df_g4_g1) <- df_g4_g1$Geneid

# Get all gene IDs in consistent order (use df_g2_g1 as reference)
all_geneids <- df_g2_g1$Geneid

# Extract counts with proper alignment
counts_mat <- data.frame(
  Geneid = all_geneids,
  X10 = df_g4_g1[all_geneids, "X10"],
  X11 = df_g4_g1[all_geneids, "X11"],
  X12 = df_g4_g1[all_geneids, "X12"],
  X7 = df_g3_g1[all_geneids, "X7"],
  X8 = df_g3_g1[all_geneids, "X8"],
  X9 = df_g3_g1[all_geneids, "X9"]
)

# For group4 vs group3 comparison, we only care about these 6 samples
group_samples <- c("X7", "X8", "X9", "X10", "X11", "X12")

# Apply filtering condition: genes with >= 5 reads in at least 3 of the 6 relevant samples
keep_condition <- rowSums(counts_mat[, group_samples] >= 5) >= 3

# Create a named logical vector for easy lookup
names(keep_condition) <- counts_mat$Geneid

cat("Total genes:", nrow(counts_mat), "\n")
cat("Genes passing filter (>=5 reads in >=3 samples):", sum(keep_condition), "\n")
cat("Genes failing filter (reads too low):", sum(!keep_condition), "\n")

# Verify the problematic gene
cat("\nVerification for ENSG00000000003:\n")
gene_test <- counts_mat[counts_mat$Geneid == "ENSG00000000003", ]
cat("X10, X11, X12:", gene_test$X10, gene_test$X11, gene_test$X12, "\n")
cat("X7, X8, X9:", gene_test$X7, gene_test$X8, gene_test$X9, "\n")

# Create metadata for the 6 samples
meta <- data.frame(
  sample_id = group_samples,
  Group = c(rep("group3", 3), rep("group4", 3))
)
rownames(meta) <- meta$sample_id

# Extract count matrix for the 6 samples only
count_matrix <- as.matrix(counts_mat[, group_samples])
rownames(count_matrix) <- counts_mat$Geneid

# Create DESeq2 dataset with ALL genes
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = meta,
  design = ~ Group
)

# Run DESeq2 analysis on ALL genes
dds <- DESeq(dds)

# Extract results for group4 vs group3 (group4 as treatment, group3 as control)
res <- results(dds, contrast = c("Group", "group4", "group3"))

# Apply lfc shrinkage using normal method
res <- lfcShrink(dds, contrast = c("Group", "group4", "group3"), type = "normal")

# Convert to dataframe - gene IDs are rownames
res_df <- as.data.frame(res)
gene_ids <- rownames(res_df)

# Add original counts back using gene_ids (with proper alignment)
count_matrix_subset <- count_matrix[gene_ids, , drop = FALSE]

# Create final result dataframe with correct column order
result_final <- data.frame(
  Geneid = gene_ids,
  X10 = count_matrix_subset[, "X10"],
  X11 = count_matrix_subset[, "X11"],  
  X12 = count_matrix_subset[, "X12"],
  X7 = count_matrix_subset[, "X7"],
  X8 = count_matrix_subset[, "X8"],
  X9 = count_matrix_subset[, "X9"],
  baseMean = res_df$baseMean,
  log2FoldChange = res_df$log2FoldChange,
  lfcSE = res_df$lfcSE,
  stat = res_df$stat,
  pvalue = res_df$pvalue,
  padj = res_df$padj
)

# Calculate FoldChange from log2FoldChange
result_final$FoldChange <- 2^result_final$log2FoldChange

# Handle infinite values and NA
result_final$FoldChange[is.infinite(result_final$FoldChange) | is.na(result_final$FoldChange)] <- ifelse(
  result_final$log2FoldChange[is.infinite(result_final$FoldChange) | is.na(result_final$FoldChange)] > 0, 
  1000, 0.001
)

# Get the keep condition for each gene in result_final
current_keep <- keep_condition[result_final$Geneid]

# Apply regulation rules from Detailed_Information.csv for genes that pass the filter
# For genes that fail the filter, mark as "reads too low to do accurate analysis"
result_final$regulation <- ifelse(
  !current_keep,
  "reads too low to do accurate analysis",
  ifelse(
    !is.na(result_final$pvalue) & result_final$pvalue < 0.05,
    ifelse(
      result_final$FoldChange >= 1.5, "Up",
      ifelse(result_final$FoldChange < (1/1.5), "Down", "No")
    ),
    "No"
  )
)

# Apply significant rules - only genes passing filter can be significant
result_final$significant <- ifelse(
  !current_keep,
  "reads too low to do accurate analysis",
  ifelse(
    !is.na(result_final$padj) & result_final$padj < 0.05 & result_final$regulation %in% c("Up", "Down"), 
    "Yes", "No"
  )
)

# Sort by padj for genes that pass filter, then by Geneid for genes that don't
# Create a sorting key: genes passing filter get their padj value, others get Inf
sort_key <- ifelse(!current_keep, Inf, ifelse(is.na(result_final$padj), Inf, result_final$padj))
result_final <- result_final[order(sort_key, result_final$Geneid), ]

# Save to CSV
write_csv(result_final, "DEG_results_group4-vs-group3.csv")

# Summary statistics
total_genes <- nrow(result_final)
passing_genes <- sum(current_keep, na.rm = TRUE)
failing_genes <- sum(!current_keep, na.rm = TRUE)
sig_genes <- sum(result_final$significant == "Yes", na.rm = TRUE)
up_genes <- sum(result_final$regulation == "Up" & result_final$significant == "Yes", na.rm = TRUE)
down_genes <- sum(result_final$regulation == "Down" & result_final$significant == "Yes", na.rm = TRUE)

cat("\nDESeq2 Analysis Summary (All Genes):\n")
cat("Total genes:", total_genes, "\n")
cat("Genes passing filter (>=5 reads in >=3 samples):", passing_genes, "\n")
cat("Genes failing filter (reads too low):", failing_genes, "\n")
cat("Significant DE genes (padj < 0.05 and Up/Down):", sig_genes, "\n")
cat("Up-regulated (FC >= 1.5 and pvalue < 0.05):", up_genes, "\n")
cat("Down-regulated (FC < 0.666... and pvalue < 0.05):", down_genes, "\n")

cat("\n✅ Analysis completed successfully!\n")
cat("Output file: DEG_results_group4-vs-group3.csv\n")