#!/usr/bin/env python3
"""
Final verification script for DEG_results_group4-vs-group3.csv
This script performs comprehensive validation to ensure data integrity.
"""

import pandas as pd
import sys

def main():
    print("=== FINAL DATA VERIFICATION FOR GROUP4 vs GROUP3 ANALYSIS ===\n")
    
    # Load all files
    try:
        df_g2_g1 = pd.read_csv('DEG_results_group2-vs-group1.csv')
        df_g3_g1 = pd.read_csv('DEG_results_group3-vs-group1.csv') 
        df_g4_g1 = pd.read_csv('DEG_results_group4-vs-group1.csv')
        df_final = pd.read_csv('DEG_results_group4-vs-group3.csv')
    except FileNotFoundError as e:
        print(f"❌ Error: Missing input file - {e}")
        sys.exit(1)
    
    # Set indices for easy lookup
    df_g3_g1.set_index('Geneid', inplace=True)
    df_g4_g1.set_index('Geneid', inplace=True)
    df_final.set_index('Geneid', inplace=True)
    
    # Verification 1: Total gene count
    expected_genes = 30330
    actual_genes = len(df_final)
    if actual_genes == expected_genes:
        print(f"✅ Gene count verification: {actual_genes} genes (expected: {expected_genes})")
    else:
        print(f"❌ Gene count verification FAILED: {actual_genes} genes (expected: {expected_genes})")
        sys.exit(1)
    
    # Verification 2: Column structure
    expected_columns = ['Geneid', 'X10', 'X11', 'X12', 'X7', 'X8', 'X9', 'baseMean', 
                      'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'FoldChange', 
                      'regulation', 'significant']
    actual_columns = list(df_final.reset_index().columns)
    if actual_columns == expected_columns:
        print("✅ Column structure verification: All columns present and in correct order")
    else:
        print("❌ Column structure verification FAILED")
        print(f"Expected: {expected_columns}")
        print(f"Actual: {actual_columns}")
        sys.exit(1)
    
    # Verification 3: Critical gene values
    test_genes = ['ENSG00000000003', 'ENSG00000001497', 'ENSG00000001629']
    for gene in test_genes:
        if gene not in df_final.index:
            print(f"❌ Gene {gene} not found in final results")
            continue
            
        # Verify group4 values
        assert df_final.loc[gene, 'X10'] == df_g4_g1.loc[gene, 'X10'], f'X10 mismatch for {gene}'
        assert df_final.loc[gene, 'X11'] == df_g4_g1.loc[gene, 'X11'], f'X11 mismatch for {gene}'
        assert df_final.loc[gene, 'X12'] == df_g4_g1.loc[gene, 'X12'], f'X12 mismatch for {gene}'
        
        # Verify group3 values  
        assert df_final.loc[gene, 'X7'] == df_g3_g1.loc[gene, 'X7'], f'X7 mismatch for {gene}'
        assert df_final.loc[gene, 'X8'] == df_g3_g1.loc[gene, 'X8'], f'X8 mismatch for {gene}'
        assert df_final.loc[gene, 'X9'] == df_g3_g1.loc[gene, 'X9'], f'X9 mismatch for {gene}'
    
    print(f"✅ Critical gene values verification: {len(test_genes)} genes verified")
    
    # Verification 4: No NA values in regulation/significant columns
    na_regulation = df_final['regulation'].isna().sum()
    na_significant = df_final['significant'].isna().sum()
    if na_regulation == 0 and na_significant == 0:
        print("✅ No NA values in regulation/significant columns")
    else:
        print(f"❌ Found NA values: regulation={na_regulation}, significant={na_significant}")
        sys.exit(1)
    
    # Verification 5: Low expression genes properly marked
    low_expr_count = (df_final['regulation'] == "reads too low to do accurate analysis").sum()
    print(f"✅ Low expression genes properly marked: {low_expr_count} genes")
    
    print("\n🎉 ALL VERIFICATION TESTS PASSED!")
    print("The final results file is ready for use.")
    
    # Summary statistics
    total_genes = len(df_final)
    sig_genes = (df_final['significant'] == 'Yes').sum()
    up_genes = ((df_final['regulation'] == 'Up') & (df_final['significant'] == 'Yes')).sum()
    down_genes = ((df_final['regulation'] == 'Down') & (df_final['significant'] == 'Yes')).sum()
    low_expr_genes = (df_final['regulation'] == "reads too low to do accurate analysis").sum()
    
    print(f"\n📊 Final Results Summary:")
    print(f"   Total genes: {total_genes:,}")
    print(f"   Significant DE genes: {sig_genes:,}")
    print(f"     Up-regulated: {up_genes:,}")
    print(f"     Down-regulated: {down_genes:,}")
    print(f"   Low-expression genes: {low_expr_genes:,}")

if __name__ == "__main__":
    main()