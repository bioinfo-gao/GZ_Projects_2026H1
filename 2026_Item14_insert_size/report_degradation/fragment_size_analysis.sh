#!/bin/bash

# Fragment size analysis for DNA degradation characterization
# This script uses fastp's overlap mechanism to estimate insert sizes
# and analyzes fragment length distributions in merged reads.

set -e

echo "Starting Fragment Size Analysis for DNA Degradation Assessment..."
echo "================================================================"

# Create results directory
RESULTS_DIR="results_fastp"
DEGRADATION_RESULTS="report_degradation/results"

mkdir -p "$DEGRADATION_RESULTS"

# Define samples
SAMPLES=("YCL_20" "YCL_21" "YCL_72" "YCL_73")
BASE_DATA_DIR="/home/gao/Dropbox/Quote_260203003/Raw_Data"

# Function to extract insert size from fastp JSON
extract_insert_size() {
    local json_file="$1"
    local sample="$2"
    
    if [ -f "$json_file" ]; then
        # Extract insert size peak from histogram
        INSERT_SIZE=$(jq '.insert_size.histogram | to_entries | map(select(.value > 100000)) | max_by(.value) | .key' "$json_file" 2>/dev/null || echo "0")
        echo "$sample: Insert size peak = $INSERT_SIZE bp"
    else
        echo "Warning: JSON file not found for $sample"
    fi
}

# Function to analyze merged read lengths
analyze_merged_lengths() {
    local merged_file="$1"
    local sample="$2"
    local output_file="$3"
    
    if [ -f "$merged_file" ]; then
        echo "Analyzing merged read lengths for $sample..."
        zcat "$merged_file" | awk 'NR%4==2 {print length}' | sort -n | uniq -c > "$output_file"
        
        # Calculate statistics
        TOTAL_READS=$(zcat "$merged_file" | awk 'NR%4==2' | wc -l)
        MEAN_LENGTH=$(zcat "$merged_file" | awk 'NR%4==2 {sum += length; count++} END {print sum/count}')
        MEDIAN_LENGTH=$(zcat "$merged_file" | awk 'NR%4==2 {print length}' | sort -n | awk 'NR==int(NR/2)')
        
        echo "$sample merged read statistics:" >> "$DEGRADATION_RESULTS/length_summary.txt"
        echo "  Total reads: $TOTAL_READS" >> "$DEGRADATION_RESULTS/length_summary.txt"
        echo "  Mean length: $MEAN_LENGTH bp" >> "$DEGRADATION_RESULTS/length_summary.txt"
        echo "  Median length: $MEDIAN_LENGTH bp" >> "$DEGRADATION_RESULTS/length_summary.txt"
        echo "" >> "$DEGRADATION_RESULTS/length_summary.txt"
    else
        echo "Warning: Merged file not found for $sample"
    fi
}

# Main analysis loop
for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample..."
    
    # Find FASTQ files
    R1_FILE="$BASE_DATA_DIR/$sample/${sample}_CKDL260002347-1A_23752VLT4_L8_1.fq.gz"
    R2_FILE="$BASE_DATA_DIR/$sample/${sample}_CKDL260002347-1A_23752VLT4_L8_2.fq.gz"
    
    if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
        echo "Error: FASTQ files not found for $sample"
        continue
    fi
    
    # Run fastp with overlap analysis (if not already done)
    if [ ! -f "$RESULTS_DIR/$sample/${sample}_fastp.json" ]; then
        mkdir -p "$RESULTS_DIR/$sample"
        echo "Running fastp for $sample..."
        fastp \
            -i "$R1_FILE" \
            -I "$R2_FILE" \
            -o /dev/null \
            -O /dev/null \
            --merge \
            --html "$RESULTS_DIR/$sample/${sample}_fastp.html" \
            --json "$RESULTS_DIR/$sample/${sample}_fastp.json" \
            --merged_out "$RESULTS_DIR/$sample/${sample}_merged.fq.gz" \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --thread 8
    else
        echo "Fastp results already exist for $sample, skipping..."
    fi
    
    # Extract insert size
    extract_insert_size "$RESULTS_DIR/$sample/${sample}_fastp.json" "$sample"
    
    # Analyze merged read lengths
    analyze_merged_lengths "$RESULTS_DIR/$sample/${sample}_merged.fq.gz" "$sample" "$DEGRADATION_RESULTS/${sample}_length_dist.txt"
    
    echo "Completed $sample"
    echo ""
done

# Generate final summary
echo "FRAGMENT SIZE ANALYSIS SUMMARY" > "$DEGRADATION_RESULTS/fragment_size_summary.txt"
echo "=============================" >> "$DEGRADATION_RESULTS/fragment_size_summary.txt"
echo "" >> "$DEGRADATION_RESULTS/fragment_size_summary.txt"

# Add insert size information from fastp logs (if available)
if [ -f "$RESULTS_DIR/YCL_20/YCL_20_fastp.json" ]; then
    echo "Insert Size Peaks (from fastp overlap analysis):" >> "$DEGRADATION_RESULTS/fragment_size_summary.txt"
    for sample in "${SAMPLES[@]}"; do
        if [ -f "$RESULTS_DIR/$sample/${sample}_fastp.json" ]; then
            INSERT_SIZE=$(jq '.insert_size.histogram | to_entries | map(select(.value > 100000)) | max_by(.value) | .key' "$RESULTS_DIR/$sample/${sample}_fastp.json" 2>/dev/null || echo "N/A")
            echo "  $sample: $INSERT_SIZE bp" >> "$DEGRADATION_RESULTS/fragment_size_summary.txt"
        fi
    done
    echo "" >> "$DEGRADATION_RESULTS/fragment_size_summary.txt"
fi

# Add merged read length statistics
if [ -f "$DEGRADATION_RESULTS/length_summary.txt" ]; then
    cat "$DEGRADATION_RESULTS/length_summary.txt" >> "$DEGRADATION_RESULTS/fragment_size_summary.txt"
fi

echo "Fragment size analysis complete!"
echo "Results available in: $DEGRADATION_RESULTS/"
echo ""
echo "Key files:"
echo "- fragment_size_summary.txt: Overall fragment size summary"
echo "- length_summary.txt: Merged read length statistics"  
echo "- {sample}_length_dist.txt: Detailed length distributions"
echo ""
echo "Interpretation guidelines:"
echo "- Consistent 58-60 bp fragments suggest systematic DNA degradation"
echo "- Variable fragment sizes might indicate random shearing or mixed contamination"
echo "- Correlation between insert size and merged read length validates analysis"