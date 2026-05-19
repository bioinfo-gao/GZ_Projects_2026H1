# Insert Size Analysis Results

## Fastp Overlap-based Analysis Results

The analysis was performed using **fastp v1.3.0** with overlap detection enabled on four sequencing samples.

### Key Results

| Sample | Insert Size Peak (bp) | Read Length | Merge Rate | Total Reads |
|--------|----------------------|-------------|------------|-------------|
| YCL_20 | 59 | 150 bp | 97.8% | 48,469,674 pairs |
| YCL_21 | 60 | 150 bp | 97.9% | 46,352,556 pairs |
| YCL_72 | 59 | 150 bp | 97.9% | 14,691,807 pairs |
| YCL_73 | 58 | 150 bp | 98.0% | 10,548,652 pairs |

### Analysis Details

- **Method**: Fastp overlap analysis (no reference genome required)
- **Read Configuration**: Paired-end 150 bp + 150 bp
- **Overlap Detection**: Enabled with 10% mismatch tolerance
- **Base Correction**: Applied in overlapped regions
- **Merge Rate**: >97% for all samples, indicating high-quality overlapping reads

### Interpretation

1. **Consistent Insert Sizes**: All four samples show very similar insert sizes (58-60 bp)
2. **Small Fragment Library**: The insert sizes are significantly smaller than the read length (150 bp), which explains the high merge rates (>97%)
3. **High Quality Data**: Excellent Q30 scores (>84% for all samples) and high merge rates indicate high-quality sequencing data
4. **Library Preparation**: The consistent small insert sizes suggest these samples were prepared as a specialized library type (possibly for specific applications like ATAC-seq or targeted sequencing)

### Technical Notes

- **Insert Size Peak**: Determined by fastp's internal algorithm based on overlapping read pairs
- **Accuracy**: For overlapping reads, insert size estimation is highly accurate
- **Limitations**: This method works best when insert size < combined read length (300 bp in this case)

### Files Generated

- **JSON Reports**: Detailed statistics in `results_fastp/{sample}/{sample}_fastp.json`
- **HTML Reports**: Interactive visualizations in `results_fastp/{sample}/{sample}_fastp.html`
- **Merged Reads**: Overlapping read pairs in `results_fastp/{sample}/{sample}_merged.fq.gz`

### Conclusion

All four samples (YCL_20, YCL_21, YCL_72, YCL_73) have remarkably consistent insert sizes around **58-60 base pairs**, indicating uniform library preparation and high-quality sequencing data suitable for downstream analysis.