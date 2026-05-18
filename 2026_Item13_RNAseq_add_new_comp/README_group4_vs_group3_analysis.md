# Group4 vs Group3 Differential Expression Analysis

## Overview
This document describes the correct procedure for performing differential expression analysis between group4 and group3 using existing DEG result files.

## Source Files Used
- `DEG_results_group2-vs-group1.csv` - Contains X4, X5, X6 (group2) and X1, X2, X3 (group1)
- `DEG_results_group3-vs-group1.csv` - Contains X7, X8, X9 (group3) and X1, X2, X3 (group1)  
- `DEG_results_group4-vs-group1.csv` - Contains X10, X11, X12 (group4) and X1, X2, X3 (group1)

## Key Requirements
1. **Filtering condition**: `rowSums(counts_mat >= 5) >= 3` (genes with ≥5 reads in at least 3 of the 6 samples)
2. **All genes retained**: Total gene count = 30,330 (same as other comparison files)
3. **Low expression marking**: Genes failing filter marked as "reads too low to do accurate analysis"
4. **Reference group**: group3 used as denominator (baseline)
5. **Significance criteria** (from Detailed_Information.csv):
   - Up-regulated: FoldChange ≥ 1.5 and pvalue < 0.05
   - Down-regulated: FoldChange < (1/1.5) and pvalue < 0.05
   - Significant: padj < 0.05 and (Up/Down)

## Critical Implementation Details

### Gene ID Alignment
The most critical step is proper gene ID alignment across files. Each file must be indexed by `Geneid` to ensure correct matching:

```r
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
```

### Verification Example
For gene ENSG00000000003:
- From `DEG_results_group4-vs-group1.csv`: X10=1604, X11=2523, X12=2374
- From `DEG_results_group3-vs-group1.csv`: X7=2493, X8=2107, X9=2076
- Final output must preserve these exact values

## MANDATORY DATA VERIFICATION STEPS

### Step 1: Pre-analysis Verification
Before running DESeq2 analysis, verify that extracted counts match source files exactly:

```bash
# Verify specific genes against source files
grep "ENSG00000000003" DEG_results_group4-vs-group1.csv | cut -d',' -f1-4
grep "ENSG00000000003" DEG_results_group3-vs-group1.csv | cut -d',' -f1-4  
grep "ENSG00000000003" DEG_results_group4-vs-group3.csv | cut -d',' -f1-7
```

Expected output verification:
- Source group4: `ENSG00000000003,1604,2523,2374,...`
- Source group3: `ENSG00000000003,2493,2107,2076,...`  
- Final output: `ENSG00000000003,1604,2523,2374,2493,2107,2076,...`

### Step 2: Post-analysis Integrity Check
After generating the final file, perform comprehensive validation:

```bash
# 1. Verify total gene count
wc -l DEG_results_group4-vs-group3.csv  # Should be 30331 (30330 genes + header)

# 2. Verify column count  
head -1 DEG_results_group4-vs-group3.csv | tr ',' '\n' | wc -l  # Should be 16 columns

# 3. Verify specific high-confidence genes
# High expression gene example
grep "ENSG00000001497" DEG_results_group4-vs-group3.csv

# Low expression gene example  
grep "ENSG00000310592" DEG_results_group4-vs-group3.csv

# 4. Verify no NA values in critical columns
cut -d',' -f15,16 DEG_results_group4-vs-group3.csv | grep -c "NA"  # Should be 0
```

### Step 3: Statistical Consistency Check
Verify that the analysis results are biologically reasonable:

```r
# Load final results
final_results <- read.csv("DEG_results_group4-vs-group3.csv")

# Check fold change distribution
summary(final_results$FoldChange[!is.infinite(final_results$FoldChange)])

# Verify significant genes have proper regulation labels
table(final_results$significant, final_results$regulation)

# Check that low-expression genes are properly marked
low_expr_genes <- final_results[final_results$regulation == "reads too low to do accurate analysis", ]
nrow(low_expr_genes)  # Should equal number of genes failing filter condition
```

### Step 4: Cross-validation with Source Data
Ensure that the original count values are preserved exactly:

```python
import pandas as pd

# Load source files
df_g3 = pd.read_csv('DEG_results_group3-vs-group1.csv', index_col='Geneid')
df_g4 = pd.read_csv('DEG_results_group4-vs-group1.csv', index_col='Geneid') 
df_final = pd.read_csv('DEG_results_group4-vs-group3.csv', index_col='Geneid')

# Test random sample of genes
test_genes = ['ENSG00000000003', 'ENSG00000001497', 'ENSG00000001629']

for gene in test_genes:
    # Verify group4 counts
    assert df_final.loc[gene, 'X10'] == df_g4.loc[gene, 'X10']
    assert df_final.loc[gene, 'X11'] == df_g4.loc[gene, 'X11'] 
    assert df_final.loc[gene, 'X12'] == df_g4.loc[gene, 'X12']
    
    # Verify group3 counts  
    assert df_final.loc[gene, 'X7'] == df_g3.loc[gene, 'X7']
    assert df_final.loc[gene, 'X8'] == df_g3.loc[gene, 'X8']
    assert df_final.loc[gene, 'X9'] == df_g3.loc[gene, 'X9']
    
print("✅ All verification tests passed!")
```

## Analysis Workflow
1. Read all three source files
2. Align gene IDs properly using Geneid as primary key
3. Extract raw counts for X7,X8,X9,X10,X11,X12
4. Apply filtering condition: `rowSums(counts >= 5) >= 3`
5. Create DESeq2 dataset with all 30,330 genes
6. Run DESeq2 analysis with group4 vs group3 contrast
7. Apply significance criteria from Detailed_Information.csv
8. Mark low-expression genes appropriately
9. Sort results by padj (significant genes first, then low-expression genes)
10. **PERFORM MANDATORY DATA VERIFICATION** (Steps 1-4 above)

## Output File Structure
- **Filename**: `DEG_results_group4-vs-group3.csv`
- **Columns**: `Geneid,X10,X11,X12,X7,X8,X9,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,FoldChange,regulation,significant`
- **Total genes**: 30,330
- **Sorting**: By padj (ascending), then by Geneid

## Common Pitfalls to Avoid
1. **Incorrect gene alignment**: Assuming row order is consistent across files (it's not!)
2. **Using wrong columns**: Ensure X7,X8,X9 come from group3 file and X10,X11,X12 from group4 file
3. **Filtering before DESeq2**: Apply filtering only for significance assignment, not for DESeq2 analysis
4. **Baseline confusion**: group3 is denominator (baseline), so FoldChange = group4/group3
5. **Skipping verification**: Never skip the mandatory data verification steps - this is critical for result integrity

## Final Results Summary
- Total genes: 30,330
- Genes passing filter (≥5 reads in ≥3 samples): 18,677
- Low-expression genes: 11,653
- Significant DE genes: 5,651
  - Up-regulated: 3,005
  - Down-regulated: 2,646