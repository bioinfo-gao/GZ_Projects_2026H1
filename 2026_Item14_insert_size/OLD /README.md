# Insert Size Analysis for Sequencing Samples

This project analyzes the insert size for four sequencing samples:
- YCL_20
- YCL_21
- YCL_72
- YCL_73

## Directory Structure

```
/home/gao/projects/2026_Item14_insert_size/
├── analyze_insert_size.sh    # Main analysis script (recommended)
├── analyze_insert_size.py    # Python alternative script
├── quick_insert_check.sh     # Quick analysis from existing BAM files
├── README.md                 # This file
└── results/                  # Output directory (created during analysis)
```

## Prerequisites

### Required Environment
The scripts require the **regular_bioinfo** mamba environment which already contains all necessary tools:
- `bwa` - for read alignment  
- `samtools` - for BAM processing
- `R` - for generating insert size distribution plots (optional)

### Activate the Environment
Before running any analysis, activate the regular_bioinfo environment:

```bash
# Using mamba/conda
mamba activate regular_bioinfo

# Or using conda directly
conda activate regular_bioinfo
```

You can verify the tools are available after activation:
```bash
which bwa      # Should show: /Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/bwa
which samtools # Should show: /Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/samtools
```

## Usage

### 1. Using the Bash Script (Recommended)

```bash
# Make sure you're in the regular_bioinfo environment
mamba activate regular_bioinfo

# Make the script executable
chmod +x /home/gao/projects/2026_Item14_insert_size/analyze_insert_size.sh

# Run the analysis (replace /path/to/reference.fasta with your reference genome)
cd /home/gao/projects/2026_Item14_insert_size
./analyze_insert_size.sh --reference /path/to/reference.fasta --threads 8
```

### 2. Using the Python Script

```bash
# Make sure you're in the regular_bioinfo environment
mamba activate regular_bioinfo

# Run the Python script (replace /path/to/reference.fasta with your reference genome)
python3 /home/gao/projects/2026_Item14_insert_size/analyze_insert_size.py \
    --reference /path/to/reference.fasta \
    --threads 8
```

### 3. Quick Check from Existing BAM Files

If you already have aligned BAM files, you can use the quick check script:

```bash
# Make sure you're in the regular_bioinfo environment
mamba activate regular_bioinfo

# Make script executable
chmod +x /home/gao/projects/2026_Item14_insert_size/quick_insert_check.sh

# Analyze existing BAM files
./quick_insert_check.sh /path/to/sample1.bam /path/to/sample2.bam

# Or if BAM files are in the results directory from previous runs
./quick_insert_check.sh
```

## Input Files

The scripts automatically detect the following FASTQ files:

- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_20/YCL_20_CKDL260002347-1A_23752VLT4_L8_1.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_20/YCL_20_CKDL260002347-1A_23752VLT4_L8_2.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_21/YCL_21_CKDL260002347-1A_23752VLT4_L8_1.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_21/YCL_21_CKDL260002347-1A_23752VLT4_L8_2.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_72/YCL_72_CKDL260002347-1A_23752VLT4_L8_1.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_72/YCL_72_CKDL260002347-1A_23752VLT4_L8_2.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_73/YCL_73_CKDL260002347-1A_23752VLT4_L8_1.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_73/YCL_73_CKDL260002347-1A_23752VLT4_L8_2.fq.gz`

## Output

Results will be saved in the `results/` subdirectory with the following files for each sample:

- `{sample}.sorted.bam` - Aligned reads in BAM format
- `{sample}.sorted.bam.bai` - BAM index file
- `{sample}_insert_stats.txt` - Basic insert size statistics
- `{sample}_insert_size.pdf` - Insert size distribution plot (if R is available)

## Notes

1. **Reference Genome**: You must provide a reference genome in FASTA format. The reference should be indexed with `bwa index` before running the analysis.

2. **Memory Requirements**: Alignment can be memory-intensive. Ensure you have sufficient RAM available.

3. **Time Estimation**: Processing time depends on the number of reads and available CPU cores.

4. **Alternative Approach**: If you already have aligned BAM files, you can skip the alignment step and directly analyze insert sizes using `samtools stats`.

5. **Environment Management**: Always ensure you are in the `regular_bioinfo` mamba environment before running any analysis scripts.

## Example Reference Genome Setup

If you don't have a pre-indexed reference genome:

```bash
# Make sure you're in the regular_bioinfo environment
mamba activate regular_bioinfo

# Download reference (example for human hg38)
wget -O hg38.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Index reference for BWA
bwa index hg38.fa

# Run analysis
./analyze_insert_size.sh --reference hg38.fa --threads 8
```

## Troubleshooting

- **"Command not found" errors**: Ensure you have activated the `regular_bioinfo` mamba environment
- **Permission denied**: Run `chmod +x analyze_insert_size.sh` to make the script executable
- **Insufficient memory**: Reduce the number of threads or use a machine with more RAM
- **Reference genome issues**: Ensure your reference genome is properly indexed with `bwa index`
- **Environment issues**: Verify your environment with `which bwa` and `which samtools`