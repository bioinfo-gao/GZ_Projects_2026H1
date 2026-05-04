#!/usr/bin/env Rscript
# env in bash 各种包装在了 DE_R45 环境 , Regular_bioinfo lacks ggrepel and ashr

# mamba activate DE_R45                 # # mamba activate regular_bioinfo
# 保留归一化、PCA、热图和 log2FC 排序，但不再输出伪造的 pvalue/padj。

# 我已经直接改了 4_run_DE_PCA_final.R (line 1)。

# 核心改动是把这份脚本改成了“无重复样本的 exploratory 模式”：
# 1 vs 1 时不再跑 dds <- DESeq(dds)，也不再用 vst(dds, blind = FALSE) 做正式差异分析，
# 而是改成 design = ~1 + estimateSizeFactors(dds)，然后用 log2(normalized count + 1) 做 PCA/heatmap，
# 并输出基于 normalized counts 的 log2FoldChange 排序结果。也就是说，你原来那三行在当前数据条件下是需要改的，尤其 DESeq(dds) 这一步不能继续保留。
# 另外也修了几个会影响你这次项目的点：
# Minxue 现在会自动生成 Minxue_T vs Minxue_C 对比，不再用脚本里原来那几个 MQ-07-99 / DMSO 的硬编码 contrasts。

# 热图部分补上了 first_comp，避免后面直接报错。
# 报告和输出文件说明也改了：无重复时不再宣称有可靠的 pvalue/padj，图也改成 exploratory fold-change 图。
# 标准化矩阵 Normalized_Counts.csv 现在会真正写出。

# 我做了两层验证：
# 语法解析已经通过。
# 在 conda run -n DE_R45 下实际执行时，脚本确实进入了 Exploratory analysis without biological replicates 分支，说明逻辑切换是对的。
# 这边执行最终停在写 PCA.pdf 时，是当前环境把目标目录挂成了只读文件系统，不是脚本统计流程本身报错。

PROJECT_DIR <- "/home/gao/projects/2026_Item11_Mingxue"
SCRIPT_DIR <- file.path(PROJECT_DIR, "scripts")

setwd(SCRIPT_DIR)
getwd()


# 跑完的输出文件：
# DEG_Test_1vsControl.csv
# DEG_Test_2vsControl.csv
# Volcano_Test_1vsControl.png
# Volcano_Test_2vsControl.png
# PCA.pdf
# Heatmap_top50_Test_1_vs_Control.pdf
# All_sample_gene_counts.tsv (拷贝的原始文件)

#!/usr/bin/env Rscript
# ==========================================================
# nf-core RNA-seq 下游分析：PCA + 差异表达 + 可视化 (针对LZJ项目修正版)
# 对比顺序: c("分组变量", "处理组/分子", "对照组/分母")
# log2FC = log2(处理组 / 对照组)
# 正数 = 在处理组中上调；负数 = 在对照组中上调
# ==========================================================

# ================= 0. 依赖加载 =================
library(DESeq2)
library(ashr)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(readr)
library(tidyr)
library(ggrepel)

# ================= 1. 路径设置 =================
META_FILE <- file.path(SCRIPT_DIR, "samples.csv")
COUNT_FILE <- file.path(PROJECT_DIR, "output_results/star_salmon/salmon.merged.gene_counts.tsv")
TPM_FILE <- file.path(PROJECT_DIR, "output_results/star_salmon/salmon.merged.gene_tpm.tsv")
OUT_DIR <- file.path(PROJECT_DIR, "Data_Analysis/DE_PCA_Results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ================= 1.5 定义数据读取目录 =================
READS_DIR <- file.path(PROJECT_DIR, "Data_Analysis/Reads")

dir.create(READS_DIR, showWarnings = FALSE, recursive = TRUE)

# ================= 2. 拷贝原始计数文件和TPM文件到目标目录 =================
# Create subdirectory for raw/normalized data

READS_DIR

# 定义需要拷贝的文件列表及其目标文件名
# [MODIFIED] Updated file list to include specific source paths
files_to_copy <- list(
  list(
    source = COUNT_FILE,
    dest = "All_sample_gene_counts.tsv",
    required = TRUE
  ),
  list(source = TPM_FILE, dest = "All_sample_gene_tpm.tsv", required = FALSE)
)

# 统一执行拷贝逻辑
for (file_info in files_to_copy) {
  src <- file_info$source
  dst_name <- file_info$dest
  is_required <- file_info$required

  if (file.exists(src)) {
    dst_path <- file.path(READS_DIR, dst_name)
    file.copy(src, dst_path, overwrite = TRUE)
    cat("✅ 文件已拷贝:", dst_name, "\n")
  } else {
    if (is_required) {
      stop("❌ 错误: 找不到必需文件: ", src)
    } else {
      cat("⚠️  警告: 可选文件不存在，已跳过: ", src, "\n")
    }
  }
}

# ================= 2.5 拷贝 QC 文件夹和生成分析报告 =================
# Define path for Gene Annotation file
# [MODIFIED] Added logic to copy gene annotation file

GENE_ANNOTATION_SRC <- "/home/gao/projects/Genes/human_Gene_annotation_20260202.xlsx"
GENE_ANNOTATION_DEST <- file.path(
  OUT_DIR,
  "human_Gene_annotation_20260202.xlsx"
)

if (file.exists(GENE_ANNOTATION_SRC)) {
  file.copy(GENE_ANNOTATION_SRC, GENE_ANNOTATION_DEST, overwrite = TRUE)
  cat("✅ 基因注释文件已拷贝:", GENE_ANNOTATION_DEST, "\n")
} else {
  cat("⚠️  警告: 基因注释文件不存在，已跳过: ", GENE_ANNOTATION_SRC, "\n")
}

# 定义 QC 文件夹源路径 (请根据实际路径调整，例如 nf-core 的 multiqc 输出或自定义 QC 文件夹)
# 假设 QC 文件夹位于项目根目录下的 output_results/multiqc 或当前目录下的 QC 文件夹
QC_SOURCE_CANDIDATES <- c(
  file.path(PROJECT_DIR, "output_results/multiqc"),
  file.path(PROJECT_DIR, "output_results/pipeline_info"),
  file.path(SCRIPT_DIR, "QC")
)

QC_SRC <- NULL
for (path in QC_SOURCE_CANDIDATES) {
  if (dir.exists(path)) {
    QC_SRC <- path
    break
  }
}

QC_SRC


if (!is.null(QC_SRC)) {
  # Define the specific destination path for QC files
  # [MODIFIED] Hardcoded destination path for QC folder
  QC_DEST_DIR <- file.path(PROJECT_DIR, "Data_Analysis/QC")

  # Create the destination directory if it doesn't exist
  dir.create(QC_DEST_DIR, showWarnings = FALSE, recursive = TRUE)

  # Copy contents of QC_SRC into QC_DEST_DIR
  # file.copy with recursive=TRUE copies the source folder INTO the destination if destination exists
  success <- file.copy(
    from = QC_SRC,
    to = QC_DEST_DIR,
    recursive = TRUE,
    overwrite = TRUE
  )

  if (any(success)) {
    cat("✅ QC 文件夹已拷贝:", QC_SRC, "->", QC_DEST_DIR, "\n")
  } else {
    cat("⚠️  警告: QC 文件夹拷贝失败\n")
  }
} else {
  cat("⚠️  警告: 未找到常见的 QC 文件夹路径，跳过拷贝\n")
}

# ================= 3. 读取并清洗元数据 =================
meta_raw <- read_csv(META_FILE)
meta_raw

# meta <- meta_raw %>%
#   select(Group, `Name in File`) 

# meta <- meta_raw %>%
#   select(Group, `Name in File`) %>%
#   rename(sample_id = `Name in File`) %>%
#   filter(!is.na(Group)) 
meta <- meta_raw %>%
  select(Group, `Name in File`) %>%
  rename(sample_id = `Name in File`) %>%
  filter(!is.na(Group)) %>%
  # ZG 本行需要按照真实组名修改 <<====================   # mutate(Group = factor(Group, levels = c("Control", "Test_1", "Test_2")))
  #mutate(Group = factor(Group, levels = c("DMSO", "MQ-07-99", "VK-8-101" , "MQ-07-81"))) 
  mutate(Group = factor(Group, levels = c("Minxue_C", "Minxue_T"))) %>%
  as.data.frame()

meta
rownames(meta) <- meta$sample_id

cat("✅ 元数据加载完成，有效样本数:", nrow(meta), "\n")

# ================= 4. 读取表达矩阵 & 预处理 =================
counts_raw <- read_tsv(COUNT_FILE, col_types = cols())
head(counts_raw)

all_sample_cols <- colnames(counts_raw)[3:ncol(counts_raw)]
all_sample_cols

# valid_samples <- all_sample_cols[!all_sample_cols %in% c("TperMix", "TtriMix")]
valid_samples <- all_sample_cols

meta <- meta[meta$sample_id %in% valid_samples, ]
meta$Group <- droplevels(meta$Group)
counts_mat <- as.matrix(counts_raw[, valid_samples])
rownames(counts_mat) <- counts_raw$gene_id
counts_mat <- counts_mat[, meta$sample_id, drop = FALSE]

counts_mat <- round(counts_mat)                                   # to integer
#keep <- rowSums(counts_mat >= 10) >= 4
keep <- rowSums(counts_mat >= 10) >= 1
  
counts_mat <- counts_mat[keep, , drop = FALSE]

head(counts_mat)
cat("✅ 表达矩阵加载完成，过滤后保留基因数:", nrow(counts_mat), "\n")

group_counts <- table(meta$Group)
has_replicates <- all(group_counts >= 2)
analysis_mode <- if (has_replicates) {
  "DESeq2 differential expression"
} else {
  "Exploratory analysis without biological replicates"
}

cat("✅ 当前分析模式:", analysis_mode, "\n")
print(group_counts)

if (!has_replicates) {
  cat(
    "⚠️  每组样本数不足 2，无法进行可靠的 DESeq2 离散度估计和正式差异检验。\n",
    "⚠️  后续仅输出归一化表达、PCA/heatmap 和基于 normalized counts 的 exploratory log2FC 排序结果。\n",
    sep = ""
  )
}



# ================= 5. 建模/归一化 =================
if (has_replicates) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = meta,
    design = ~Group
  )
  dds <- DESeq(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  vsd <- vst(dds, blind = FALSE)
  expr_mat <- assay(vsd)
} else {
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = meta,
    design = ~1
  )
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  expr_mat <- log2(norm_counts + 1)
}


# ================= 6. PCA 分析 =================
n_top_var_genes <- min(500, nrow(expr_mat))
rv <- apply(expr_mat, 1, var)
select_genes <- order(rv, decreasing = TRUE)[seq_len(n_top_var_genes)]
pca <- prcomp(t(expr_mat[select_genes, , drop = FALSE]), scale. = FALSE)
percentVarAll <- 100 * (pca$sdev ^ 2 / sum(pca$sdev ^ 2))
percentVar <- c(round(percentVarAll[1], 1), 0)
if (length(percentVarAll) >= 2) {
  percentVar[2] <- round(percentVarAll[2], 1)
}

pca_data <- data.frame(
  name = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = if (ncol(pca$x) >= 2) pca$x[, 2] else 0,
  Group = meta$Group[match(rownames(pca$x), meta$sample_id)]
)

if (any(is.na(pca_data$Group))) {
  stop("❌ 错误: PCA 样本名与 meta$sample_id 未正确匹配，无法作图。")
}

pca_data 
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = Group, label = name)) +
  geom_point(size = 2.5, alpha = 0.9) + # 🔻 增大点大小和透明度以提高可见性
  geom_text_repel(
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = 20
  ) + # 🔻 使用 geom_text_repel 避免标签与点重合
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 12) + # 🔻 恢复全局字体大小以提高可读性
  scale_color_brewer(palette = "Set1") + # 🔻 使用对比度更高的颜色 palette
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )
print(OUT_DIR)

ggsave(file.path(OUT_DIR, "PCA.pdf"), p_pca, width = 8, height = 6, dpi = 300)
cat("✅ PCA plot saved\n")

# ================= 7. 差异表达/探索性比较 =================
# ✅ 对比组定义：c("分组列名", "处理组/分子", "对照组/分母")
# log2FC = log2(处理组均值 / 对照组均值)
# 正数 = 在处理组(第一个)中上调；负数 = 在对照组(第二个)中上调

group_levels_in_use <- levels(droplevels(meta$Group))

if (length(group_levels_in_use) < 2) {
  stop("❌ 错误: 至少需要 2 个分组才能进行比较。")
}

if (length(group_levels_in_use) == 2) {
  contrasts <- list(c("Group", group_levels_in_use[2], group_levels_in_use[1]))
} else {
  baseline_group <- group_levels_in_use[1]
  contrasts <- lapply(group_levels_in_use[-1], function(grp) {
    c("Group", grp, baseline_group)
  })
}

res_list <- list()
sig_col_name <- if (has_replicates) {
  "sig (padj<0.05 & |log2FC|>=1)"
} else {
  "sig (|log2FC|>=1)"
}

for (comp in contrasts) {
  grp_treatment <- comp[2]
  grp_control <- comp[3]
  comp_name <- paste(grp_treatment, "vs", grp_control, sep = "_")
  cat(paste0(
    "\n🔍 正在分析: ",
    comp_name,
    " (处理:",
    grp_treatment,
    " vs 对照:",
    grp_control,
    ") [",
    analysis_mode,
    "]\n"
  ))

  samples_treatment <- meta$sample_id[meta$Group == grp_treatment]
  samples_control <- meta$sample_id[meta$Group == grp_control]
  raw_sub <- counts_raw[,
    c("gene_id", "gene_name", samples_treatment, samples_control),
    drop = FALSE
  ]

  if (has_replicates) {
    res <- lfcShrink(dds, contrast = comp, type = "ashr")
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df$mean_treatment_norm <- rowMeans(
      norm_counts[, samples_treatment, drop = FALSE]
    )
    res_df$mean_control_norm <- rowMeans(
      norm_counts[, samples_control, drop = FALSE]
    )
  } else {
    mean_treatment_norm <- rowMeans(
      norm_counts[, samples_treatment, drop = FALSE]
    )
    mean_control_norm <- rowMeans(
      norm_counts[, samples_control, drop = FALSE]
    )

    res_df <- data.frame(
      gene_id = rownames(norm_counts),
      baseMean = rowMeans(norm_counts),
      log2FoldChange = log2((mean_treatment_norm + 0.5) / (mean_control_norm + 0.5)),
      lfcSE = NA_real_,
      pvalue = NA_real_,
      padj = NA_real_,
      mean_treatment_norm = mean_treatment_norm,
      mean_control_norm = mean_control_norm,
      stringsAsFactors = FALSE
    )
  }

  res_df <- left_join(res_df, raw_sub, by = "gene_id")

  final_cols <- c(
    "gene_id",
    "gene_name",
    samples_treatment,
    samples_control,
    "mean_treatment_norm",
    "mean_control_norm",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "pvalue",
    "padj"
  )
  res_df <- res_df[, final_cols]

  if (has_replicates) {
    res_df <- res_df %>% arrange(padj, desc(abs(log2FoldChange)))
    res_df[[sig_col_name]] <- case_when(
      !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 1 ~ "Up",
      !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -1 ~ "Down",
      TRUE ~ "NS"
    )
  } else {
    res_df <- res_df %>%
      mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
      arrange(desc(abs_log2FoldChange), desc(baseMean))
    res_df[[sig_col_name]] <- case_when(
      res_df$log2FoldChange >= 1 ~ "Up",
      res_df$log2FoldChange <= -1 ~ "Down",
      TRUE ~ "NS"
    )
  }

  write_csv(res_df, file.path(OUT_DIR, paste0("DEG_", comp_name, ".csv")))
  res_list[[comp_name]] <- res_df

  if (has_replicates) {
    res_df$negLog10Padj <- -log10(res_df$padj)
    res_df$negLog10Padj[!is.finite(res_df$negLog10Padj)] <- NA

    top_labels <- res_df %>%
      filter(.data[[sig_col_name]] != "NS", !is.na(negLog10Padj)) %>%
      arrange(padj, desc(abs(log2FoldChange))) %>%
      head(10) %>%
      mutate(
        label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)
      )

    p_vol <- ggplot(
      res_df,
      aes(x = log2FoldChange, y = negLog10Padj, color = .data[[sig_col_name]])
    ) +
      geom_point(alpha = 0.7, size = 0.5) +
      scale_color_manual(
        values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey80"),
        labels = c(
          "Up" = "Upregulated",
          "Down" = "Downregulated",
          "NS" = "Not Significant"
        )
      ) +
      theme_bw(base_size = 10) +
      labs(
        title = paste(
          grp_treatment,
          "vs",
          grp_control,
          "(log2FC > 0 =",
          grp_treatment,
          "upregulated)"
        ),
        x = "log2 Fold Change",
        y = "-log10(adj. P-value)"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9)
      ) +
      geom_text_repel(
        data = top_labels,
        aes(label = label),
        size = 2,
        box.padding = 0.3,
        max.overlaps = 20,
        color = "black",
        fontface = "plain"
      )

    ggsave(
      file.path(OUT_DIR, paste0("Volcano_", comp_name, ".png")),
      p_vol,
      width = 8,
      height = 6,
      dpi = 300
    )
  } else {
    res_df$log10BaseMean <- log10(res_df$baseMean + 1)

    top_labels <- res_df %>%
      filter(.data[[sig_col_name]] != "NS") %>%
      arrange(desc(abs(log2FoldChange)), desc(baseMean)) %>%
      head(10) %>%
      mutate(
        label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)
      )

    p_fc <- ggplot(
      res_df,
      aes(x = log2FoldChange, y = log10BaseMean, color = .data[[sig_col_name]])
    ) +
      geom_point(alpha = 0.7, size = 0.5) +
      scale_color_manual(
        values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey80"),
        labels = c(
          "Up" = "Higher in treatment",
          "Down" = "Higher in control",
          "NS" = "Not Highlighted"
        )
      ) +
      theme_bw(base_size = 10) +
      labs(
        title = paste(
          grp_treatment,
          "vs",
          grp_control,
          "(exploratory; no replicates)"
        ),
        x = "log2 Fold Change",
        y = "log10(mean normalized count + 1)"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9)
      ) +
      geom_text_repel(
        data = top_labels,
        aes(label = label),
        size = 2,
        box.padding = 0.3,
        max.overlaps = 20,
        color = "black",
        fontface = "plain"
      )

    ggsave(
      file.path(OUT_DIR, paste0("Exploratory_", comp_name, ".png")),
      p_fc,
      width = 8,
      height = 6,
      dpi = 300
    )
  }

  cat(paste0(
    "✅ ",
    comp_name,
    " completed, highlighted genes: ",
    sum(res_df[[sig_col_name]] != "NS"),
    "\n"
  ))
}

# ================= 8. 热图 =================
first_comp <- names(res_list)[1]
top_genes_df <- res_list[[first_comp]] %>%
  filter(.data[[sig_col_name]] != "NS")

if (nrow(top_genes_df) == 0) {
  top_genes_df <- res_list[[first_comp]] %>%
    arrange(desc(abs(log2FoldChange)), desc(baseMean)) %>%
    head(50)
} else if (has_replicates) {
  top_genes_df <- top_genes_df %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(50)
} else {
  top_genes_df <- top_genes_df %>%
    arrange(desc(abs(log2FoldChange)), desc(baseMean)) %>%
    head(50)
}

top_genes_ids <- top_genes_df$gene_id

mat <- expr_mat[top_genes_ids, , drop = FALSE]
mat <- t(scale(t(mat)))
mat[is.na(mat)] <- 0

# ✅ 将行名替换为 gene_name，如果 gene_name 为空则使用 gene_id
gene_names_for_plot <- top_genes_df$gene_name
# 处理缺失或空的 gene_name
gene_names_for_plot[
  is.na(gene_names_for_plot) | gene_names_for_plot == ""
] <- top_genes_df$gene_id[
  is.na(gene_names_for_plot) | gene_names_for_plot == ""
]
# 确保行名唯一，如果有重复，可以添加后缀或保留原ID，这里简单处理：
# 如果存在重复的 gene_name，pheatmap 可能会报错或显示不全。
# 为了安全，我们检查重复并必要时回退到 ID
if (any(duplicated(gene_names_for_plot))) {
  # 简单的去重策略：如果重复，追加 gene_id 以确保唯一性
  gene_names_for_plot[duplicated(gene_names_for_plot)] <- paste0(
    gene_names_for_plot[duplicated(gene_names_for_plot)],
    "_",
    top_genes_df$gene_id[duplicated(gene_names_for_plot)]
  )
}

rownames(mat) <- gene_names_for_plot

# 创建正确的annotation数据框，确保是字符类型而不是因子
annotation_df <- data.frame(
  Group = as.character(meta$Group),
  row.names = meta$sample_id
)

pheatmap(
  mat,
  annotation_col = annotation_df,
  filename = file.path(OUT_DIR, "heatmap.pdf"),
  show_rownames = TRUE,
  main = paste("Top 50 genes:", first_comp),
  fontsize = 7,
  fontsize_row = 5,
  fontsize_col = 7,
  fontfamily = "sans",
  legend_labels = c("Low", "High")
) # Ensure legend labels are English if auto-generated

cat("✅ Heatmap generated\n")

# ================= 9. 保存标准化计数矩阵 =================
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$gene_id <- rownames(norm_counts_df)
gene_names <- counts_raw[, c("gene_id", "gene_name")]
norm_counts_df <- left_join(norm_counts_df, gene_names, by = "gene_id")
final_norm_cols <- c("gene_id", "gene_name", colnames(norm_counts))
norm_counts_df <- norm_counts_df[, final_norm_cols]
write_csv(norm_counts_df, file.path(OUT_DIR, "Normalized_Counts.csv"))
cat("✅ 标准化计数矩阵已保存\n")


# ================= 10. 生成分析报告 =================
# Generate Bioinformatics_Analysis_Report.md with detailed content
report_file <- file.path(OUT_DIR, "Bioinformatics_Analysis_Report.md")


# Calculate summary statistics for the report
deg_summary <- lapply(res_list, function(df) {
  up <- sum(df[[sig_col_name]] == "Up", na.rm = TRUE)
  down <- sum(df[[sig_col_name]] == "Down", na.rm = TRUE)
  total_sig <- up + down
  list(up = up, down = down, total = total_sig)
})

# Build report content
report_content <- c(
  "# Bioinformatics Analysis Report",
  "",
  paste("Date:", Sys.Date()),
  paste("Project:", "2026_Item11_Mingxue"),
  "",
  "## 1. Overview",
  "This report summarizes the differential expression analysis and quality control metrics for the RNA-seq dataset.",
  paste("- **Analysis Mode**:", analysis_mode),
  ifelse(
    has_replicates,
    "- **Analysis Tool**: DESeq2 with shrinkage of log2 fold changes",
    "- **Analysis Tool**: DESeq2 size-factor normalization only; no formal hypothesis testing"
  ),
  ifelse(
    has_replicates,
    "- **Normalization**: VST (Variance Stabilizing Transformation) for PCA/Heatmap, Median-of-ratios for DE",
    "- **Normalization**: Median-of-ratios normalization, followed by log2(normalized count + 1) for PCA/Heatmap"
  ),
  ifelse(
    has_replicates,
    "- **Significance Thresholds**: padj < 0.05, |log2FoldChange| >= 1",
    "- **Highlight Thresholds**: |log2FoldChange| >= 1 (exploratory only; p-values/padj are not reported)"
  ),
  "",
  "## 2. Quality Control (QC)",
  ifelse(
    !is.null(QC_SRC),
    c(
      "- QC reports were generated using MultiQC.",
      "- Raw data quality and alignment metrics are available in the `QC/` directory.",
      "- Please refer to `QC/multiqc_report.html` for detailed interactive plots."
    ),
    c(
      "- ⚠️ QC folder not found or not copied.",
      "- Ensure raw data QC was performed prior to this step."
    )
  ),
  "",
  "## 3. Differential Expression Analysis Results",
  ""
)

# Add DEG statistics for each contrast
for (name in names(deg_summary)) {
  stats <- deg_summary[[name]]
  report_content <- c(
    report_content,
    paste0("### Contrast: ", name),
    paste0("- Total Highlighted Genes: ", stats$total),
    paste0("  - Upregulated: ", stats$up),
    paste0("  - Downregulated: ", stats$down),
    paste0("- Output File: `DEG_", name, ".csv`"),
    ""
  )
}

report_content <- c(
  report_content,
  "## 4. Visualizations",
  "",
  "### Principal Component Analysis (PCA)",
  "- **File**: `PCA.pdf`",
  "- **Description**: Shows sample clustering based on the top variable genes used for PCA.",
  "",
  "### Fold-change Plot",
  ifelse(
    has_replicates,
    "- **Files**: `Volcano_*.png`",
    "- **Files**: `Exploratory_*.png`"
  ),
  ifelse(
    has_replicates,
    "- **Description**: Displays the relationship between statistical significance (-log10 padj) and magnitude of change (log2FC).",
    "- **Description**: Exploratory scatter plot showing log2FC against mean normalized expression."
  ),
  "",
  "### Heatmap",
  "- **File**: `heatmap.pdf`",
  ifelse(
    has_replicates,
    "- **Description**: Hierarchical clustering of the top 50 differentially expressed genes across all samples.",
    "- **Description**: Hierarchical clustering of the top 50 highest-change genes across all samples."
  ),
  "",
  "## 5. Generated Data Files",
  "",
  "| File Name | Description |",
  "| :--- | :--- |",
  "| `All_sample_gene_counts.tsv` | Raw count matrix for all samples. |",
  "| `All_sample_gene_tpm.tsv` | TPM matrix, if available. |",
  "| `Normalized_Counts.csv` | DESeq2 normalized counts for downstream analysis. |",
  ifelse(
    has_replicates,
    "| `DEG_*.csv` | Differential expression results including log2FC, p-values, and base means. |",
    "| `DEG_*.csv` | Exploratory fold-change table including normalized means, log2FC, and base means. |"
  ),
  "| `PCA.pdf` | PCA plot showing sample relationships. |",
  ifelse(
    has_replicates,
    "| `Volcano_*.png` | Volcano plots for each contrast. |",
    "| `Exploratory_*.png` | Exploratory fold-change scatter plots for each contrast. |"
  ),
  "| `heatmap.pdf` | Heatmap of top genes. |",
  ifelse(
    !is.null(QC_SRC),
    "| `QC/` | Directory containing MultiQC and other QC reports. |",
    "| `QC/` | Not available. |"
  )
)

writeLines(report_content, con = report_file)
cat("✅ 详细分析报告已生成:", report_file, "\n")

cat("\n🎉 全部分析完成！结果已保存至:", OUT_DIR, "\n")
