# Insert Size Analysis for Sequencing Samples

This project provides two methods to analyze insert sizes for four sequencing samples:
1. **Fastp-based Overlap Analysis**: A rapid, reference-free method suitable for quality control and quick estimation.
2. **Reference-based Alignment Analysis**: A precise method using BWA and Samtools, requiring a reference genome.

## Target Samples

- YCL_20
- YCL_21  
- YCL_72
- YCL_73

## Directory Structure

```
/home/gao/projects/2026_Item14_insert_size/
├── analyze_with_fastp.sh       # Fastp-based analysis script (Recommended for quick check)
├── analyze_insert_size.sh      # Reference-based analysis script
├── analyze_insert_size.py      # Python alternative for reference-based analysis
├── quick_insert_check.sh       # Quick analysis from existing BAM files
├── parse_fastp_results.py      # Utility to parse fastp JSON results
├── README.md                   # This file
├── results/                    # Output directory for reference-based analysis
└── fastp_results/              # Output directory for fastp-based analysis
```

## Method 1: Fastp-based Insert Size Analysis (Recommended for Quick Check)

This method uses **fastp's overlap mechanism** to estimate insert sizes without requiring reference genome alignment. It is significantly faster and ideal for initial quality control.

### How Fastp Overlap Analysis Works

Fastp estimates insert sizes by analyzing overlapping regions between paired-end reads:

- **When insert size < read length × 2**: Reads overlap and fastp can directly measure the overlap length.
- **Insert size estimation**: `insert_size ≈ read1_length + read2_length - overlap_length`
- **When no overlap detected**: Insert size is larger than combined read length (provides lower bound only).

### Environment Requirements

The analysis requires the **regular_bioinfo** mamba environment which contains `fastp`:

```bash
# Activate the environment
mamba activate regular_bioinfo

# Verify fastp is available
which fastp
```

### Usage

#### 1. Run the Analysis Script

```bash
# Activate environment
mamba activate regular_bioinfo

# Make script executable
chmod +x /home/gao/projects/2026_Item14_insert_size/analyze_with_fastp.sh

# Run analysis
cd /home/gao/projects/2026_Item14_insert_size
./analyze_with_fastp.sh
```

#### 2. View Results

Results will be saved in `/home/gao/projects/2026_Item14_insert_size/fastp_results/`:

- `{sample}_fastp.json` - Detailed JSON report with overlap statistics
- `{sample}_fastp.html` - Interactive HTML report  
- `summary.txt` - Consolidated summary of all samples

### Output Interpretation

#### Key Metrics in Summary

- **Total Reads**: Number of read pairs processed
- **Overlapped Reads**: Number of read pairs with detectable overlap
- **Overlap Rate**: Percentage of reads with overlap
- **Estimated Insert Size**: Calculated insert size for overlapped reads
- **Average Overlap Length**: Mean overlap length in base pairs

#### Limitations

- **Only accurate for small insert sizes**: When insert size < combined read length.
- **Lower bound for large inserts**: When no overlap detected, only provides minimum possible insert size.
- **Read length dependent**: Accuracy depends on actual read lengths.

### Alternative: Parse Results with Python

For programmatic analysis of results:

```bash
# Parse JSON reports into summary table
python3 parse_fastp_results.py --input-dir fastp_results --output insert_size_summary.txt
```

---

## Method 2: Reference-based Alignment Analysis (Precise)

This method aligns reads to a reference genome using BWA and calculates insert sizes from the resulting BAM files using Samtools. It provides precise insert size distributions for all fragments, regardless of size.

### Prerequisites

#### Required Environment
The scripts require the **regular_bioinfo** mamba environment which already contains all necessary tools:
- `bwa` - for read alignment  
- `samtools` - for BAM processing
- `R` - for generating insert size distribution plots (optional)

#### Activate the Environment
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

### Usage

#### 1. Using the Bash Script

```bash
# Make sure you're in the regular_bioinfo environment
mamba activate regular_bioinfo

# Make the script executable
chmod +x /home/gao/projects/2026_Item14_insert_size/analyze_insert_size.sh

# Run the analysis (replace /path/to/reference.fasta with your reference genome)
cd /home/gao/projects/2026_Item14_insert_size
./analyze_insert_size.sh --reference /path/to/reference.fasta --threads 8
```

#### 2. Using the Python Script

```bash
# Make sure you're in the regular_bioinfo environment
mamba activate regular_bioinfo

# Run the Python script (replace /path/to/reference.fasta with your reference genome)
python3 /home/gao/projects/2026_Item14_insert_size/analyze_insert_size.py \
    --reference /path/to/reference.fasta \
    --threads 8
```

#### 3. Quick Check from Existing BAM Files

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

### Input Files

The scripts automatically detect the following FASTQ files:

- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_20/YCL_20_CKDL260002347-1A_23752VLT4_L8_1.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_20/YCL_20_CKDL260002347-1A_23752VLT4_L8_2.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_21/YCL_21_CKDL260002347-1A_23752VLT4_L8_1.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_21/YCL_21_CKDL260002347-1A_23752VLT4_L8_2.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_72/YCL_72_CKDL260002347-1A_23752VLT4_L8_1.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_72/YCL_72_CKDL260002347-1A_23752VLT4_L8_2.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_73/YCL_73_CKDL260002347-1A_23752VLT4_L8_1.fq.gz`
- `/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_73/YCL_73_CKDL260002347-1A_23752VLT4_L8_2.fq.gz`

### Output

Results will be saved in the `results/` subdirectory with the following files for each sample:

- `{sample}.sorted.bam` - Aligned reads in BAM format
- `{sample}.sorted.bam.bai` - BAM index file
- `{sample}_insert_stats.txt` - Basic insert size statistics
- `{sample}_insert_size.pdf` - Insert size distribution plot (if R is available)

### Notes

1. **Reference Genome**: You must provide a reference genome in FASTA format. The reference should be indexed with `bwa index` before running the analysis.

2. **Memory Requirements**: Alignment can be memory-intensive. Ensure you have sufficient RAM available.

3. **Time Estimation**: Processing time depends on the number of reads and available CPU cores.

4. **Alternative Approach**: If you already have aligned BAM files, you can skip the alignment step and directly analyze insert sizes using `samtools stats`.

5. **Environment Management**: Always ensure you are in the `regular_bioinfo` mamba environment before running any analysis scripts.

### Example Reference Genome Setup

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
- **Permission denied**: Run `chmod +x analyze_insert_size.sh` or `chmod +x analyze_with_fastp.sh` to make the script executable
- **Insufficient memory**: Reduce the number of threads or use a machine with more RAM
- **Reference genome issues**: Ensure your reference genome is properly indexed with `bwa index`
- **Environment issues**: Verify your environment with `which bwa`, `which samtools`, and `which fastp`