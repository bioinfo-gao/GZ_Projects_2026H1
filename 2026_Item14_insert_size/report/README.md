# Insert Size Analysis Report Package

This directory contains all key code and evidence for the primer dimer contamination analysis of four sequencing samples (YCL_20, YCL_21, YCL_72, YCL_73).

## File Structure

```
report/
├── FINAL_ANALYSIS_REPORT.md          # Comprehensive analysis report with executive summary
├── README.md                         # This file
├── analyze_with_fastp.sh             # Main fastp-based insert size analysis pipeline
├── analyze_gc_content.py             # GC content validation script
├── analyze_sequence_complexity.py    # Detailed sequence complexity analysis
└── analyze_first_50bp.py             # Specialized analysis of first 50bp sequences
```

## Key Evidence Chain

The analysis demonstrates **primer dimer contamination** through multiple independent lines of evidence:

1. **Physical Size Violation**: Insert sizes of 58-60 bp vs expected 500 bp
2. **GC Content Abnormality**: 60-65% GC vs expected 40-45% for mammalian DNA  
3. **Sequence Complexity**: Extremely low k-mer diversity (0.1%) despite high sequence uniqueness
4. **Adapter Contamination**: Direct detection of Illumina adapter sequences

## How to Reproduce Analysis

### Prerequisites
- Conda environment: `regular_bioinfo` with fastp installed
- Input data in: `/home/gao/Dropbox/Quote_260203003/Raw_Data/`

### Steps
1. **Activate environment**:
   ```bash
   mamba activate regular_bioinfo
   ```

2. **Run main analysis**:
   ```bash
   cd /home/gao/projects/2026_Item14_insert_size/report
   chmod +x analyze_with_fastp.sh
   ./analyze_with_fastp.sh
   ```

3. **Validate GC content**:
   ```bash
   python3 analyze_gc_content.py --results-dir ../results_fastp
   ```

4. **Analyze sequence complexity**:
   ```bash
   python3 analyze_sequence_complexity.py --results-dir ../results_fastp
   python3 analyze_first_50bp.py --results-dir ../results_fastp
   ```

## Critical Findings Summary

- **All four samples** show identical contamination patterns
- **Data is unsuitable** for intended biological analysis  
- **Immediate action required**: Stop downstream analysis and re-prepare libraries
- **Root cause**: Likely size selection failure during library preparation

## For Management Reporting

The `FINAL_ANALYSIS_REPORT.md` file contains:
- Executive summary suitable for management review
- Technical details for bioinformatics team
- Specific recommendations for next steps
- Complete evidence chain with statistical support

This package provides complete transparency and reproducibility for the quality control findings.