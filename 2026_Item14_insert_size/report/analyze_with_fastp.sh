#!/bin/bash

# Fastp-based Insert Size Analysis for Primer Dimer Detection
# This script uses fastp's overlap mechanism to estimate insert sizes
# without requiring a reference genome.

set -e

# Activate the regular_bioinfo environment
echo "Activating regular_bioinfo environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate regular_bioinfo

# Verify fastp is available
if ! command -v fastp &> /dev/null; then
    echo "Error: fastp not found in PATH"
    echo "Please ensure you are in the regular_bioinfo environment"
    exit 1
fi

echo "Starting Fastp Insert Size Analysis..."
echo "======================================="

# Create results directory
OUTPUT_DIR="/home/gao/projects/2026_Item14_insert_size/results_fastp"
mkdir -p "$OUTPUT_DIR"

# Define samples
SAMPLES=("YCL_20" "YCL_21" "YCL_72" "YCL_73")
BASE_DATA_DIR="/home/gao/Dropbox/Quote_260203003/Raw_Data"

# Process each sample
for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample..."
    
    sample_output_dir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_output_dir"
    
    # Find FASTQ files
    # Using wildcard to find R1 and R2 files, assuming standard naming conventions within sample folders
    R1_FILE=$(find "$BASE_DATA_DIR/$sample" -name "*_1.fq.gz" -type f | head -n 1)
    R2_FILE=$(find "$BASE_DATA_DIR/$sample" -name "*_2.fq.gz" -type f | head -n 1)
    
    if [ -z "$R1_FILE" ] || [ -z "$R2_FILE" ] || [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
        echo "Error: FASTQ files not found for $sample in $BASE_DATA_DIR/$sample"
        exit 1
    fi
    
    echo "  R1: $(basename "$R1_FILE")"
    echo "  R2: $(basename "$R2_FILE")"
    
    # Run fastp with overlap analysis
    # --disable_adapter_trimming, --disable_quality_filtering, --disable_length_filtering:
    # Disable filtering to keep all reads for accurate insert size distribution estimation.
    # --merge: Enable merging to detect overlaps.
    # -o /dev/null -O /dev/null: Discard unmerged reads output to save space, as we focus on insert size stats.
    
    echo "  Running fastp overlap analysis..."
    fastp \
        -i "$R1_FILE" \
        -I "$R2_FILE" \
        -o /dev/null \
        -O /dev/null \
        --merge \
        --html "$sample_output_dir/${sample}_fastp.html" \
        --json "$sample_output_dir/${sample}_fastp.json" \
        --merged_out "$sample_output_dir/${sample}_merged.fq.gz" \
        --disable_adapter_trimming \
        --disable_quality_filtering \
        --disable_length_filtering \
        --thread 8 \
        --verbose
    
    echo "Completed $sample"
    echo ""
done

# Generate summary report
echo "Generating summary report..."
generate_summary() {
    SUMMARY_FILE="$OUTPUT_DIR/summary.txt"
    echo "Fastp Insert Size Analysis Summary" > "$SUMMARY_FILE"
    echo "==================================" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"
    
    for sample in "${SAMPLES[@]}"; do
        json_file="$OUTPUT_DIR/$sample/${sample}_fastp.json"
        
        if [ -f "$json_file" ]; then
            echo "Sample: $sample" >> "$SUMMARY_FILE"
            
            # Extract key metrics using jq
            total_reads=$(jq -r '.summary.before_filtering.total_reads // 0' "$json_file")
            merged_reads=$(jq -r '.merge.merged_reads // 0' "$json_file")
            merge_rate=$(jq -r '.merge.merge_rate // 0' "$json_file")
            
            echo "  Total reads: $((total_reads / 2))" >> "$SUMMARY_FILE"
            echo "  Merged reads: $merged_reads" >> "$SUMMARY_FILE"
            echo "  Merge rate: $(printf "%.2f%%" $(echo "$merge_rate * 100" | bc -l))" >> "$SUMMARY_FILE"
            
            # Estimate insert size from read lengths
            r1_len=$(jq -r '.summary.before_filtering.read1_mean_length // 0' "$json_file")
            r2_len=$(jq -r '.summary.before_filtering.read2_mean_length // 0' "$json_file")
            
            if [ "$merged_reads" -gt 0 ]; then
                # For merged reads, insert size ≈ read1_len + read2_len - overlap_len
                # average_merged_length represents the length of the merged read, which is effectively the insert size for fully merged reads
                # However, strictly speaking: Insert Size = R1 Len + R2 Len - Overlap Len.
                # Merged Read Length = R1 Len + R2 Len - Overlap Len.
                # So Merged Read Length IS the Insert Size for those merged reads.
                avg_merged_len=$(jq -r '.merge.average_merged_length // 0' "$json_file")
                if [ "$avg_merged_len" -gt 0 ]; then
                    echo "  Estimated insert size (from merged reads): ~$avg_merged_len bp" >> "$SUMMARY_FILE"
                fi
            else
                min_insert=$((r1_len + r2_len))
                echo "  No merged reads detected" >> "$SUMMARY_FILE"
                echo "  Minimum insert size: > $min_insert bp" >> "$SUMMARY_FILE"
            fi
            
            echo "" >> "$SUMMARY_FILE"
        else
            echo "  WARNING: JSON report not found for $sample" >> "$SUMMARY_FILE"
            echo "" >> "$SUMMARY_FILE"
        fi
    done
    
    echo "Summary report generated: $SUMMARY_FILE"
}

# Check if jq is available for JSON parsing
if command -v jq &> /dev/null; then
    generate_summary
else
    echo "Warning: jq not available, skipping detailed summary generation"
    echo "Install jq to get detailed JSON parsing: conda install -c conda-forge jq"
fi

echo ""
echo "Analysis complete! Results in: $OUTPUT_DIR/"
echo ""
echo "Key metrics to check:"
echo "- Insert size peak in JSON files (merge.average_merged_length)"
echo "- Merge rate in HTML reports" 
echo "- Merged read lengths and GC content"