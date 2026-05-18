#!/bin/bash

# Quick Insert Size Check Script
# This script provides multiple ways to check insert size:
# 1. From existing BAM files (fastest)
# 2. Using samtools stats on BAM files
# 3. Basic statistics from alignment data

set -e

OUTPUT_DIR="/home/gao/projects/2026_Item14_insert_size/results"
mkdir -p "$OUTPUT_DIR"

# Function to analyze insert size from BAM file using samtools stats
analyze_bam_insert_size() {
    local bam_file="$1"
    local sample_name="$2"
    
    if [ ! -f "$bam_file" ]; then
        echo "BAM file not found: $bam_file"
        return 1
    fi
    
    echo "Analyzing insert size for $sample_name from $bam_file"
    
    # Use samtools stats to get insert size information
    samtools stats "$bam_file" > "$OUTPUT_DIR/${sample_name}_samtools_stats.txt"
    
    # Extract insert size summary
    echo "=== Insert Size Summary for $sample_name ===" > "$OUTPUT_DIR/${sample_name}_insert_summary.txt"
    echo "Total reads processed:" >> "$OUTPUT_DIR/${sample_name}_insert_summary.txt"
    grep "^SN" "$OUTPUT_DIR/${sample_name}_samtools_stats.txt" | grep -E "(reads|pairs)" >> "$OUTPUT_DIR/${sample_name}_insert_summary.txt"
    
    echo "" >> "$OUTPUT_DIR/${sample_name}_insert_summary.txt"
    echo "Insert size distribution (top 10):" >> "$OUTPUT_DIR/${sample_name}_insert_summary.txt"
    grep "^IS" "$OUTPUT_DIR/${sample_name}_samtools_stats.txt" | head -10 >> "$OUTPUT_DIR/${sample_name}_insert_summary.txt"
    
    # Calculate mean and median insert size
    echo "" >> "$OUTPUT_DIR/${sample_name}_insert_summary.txt"
    echo "Calculated statistics:" >> "$OUTPUT_DIR/${sample_name}_insert_summary.txt"
    
    # Extract insert size data and calculate statistics using awk
    awk '
    /^IS/ {
        insert_size = $2
        count = $3
        total += insert_size * count
        total_reads += count
        for(i=1; i<=count; i++) {
            sizes[++size_index] = insert_size
        }
    }
    END {
        if(total_reads > 0) {
            mean = total / total_reads
            # Sort array for median calculation
            n = asort(sizes)
            if(n % 2 == 1) {
                median = sizes[int(n/2) + 1]
            } else {
                median = (sizes[n/2] + sizes[n/2 + 1]) / 2
            }
            printf "Mean insert size: %.2f\n", mean
            printf "Median insert size: %.2f\n", median
            printf "Total paired reads: %d\n", total_reads
        } else {
            print "No insert size data found"
        }
    }' "$OUTPUT_DIR/${sample_name}_samtools_stats.txt" >> "$OUTPUT_DIR/${sample_name}_insert_summary.txt"
    
    echo "Results saved to: $OUTPUT_DIR/${sample_name}_insert_summary.txt"
}

# Function to check if BAM files exist for samples
check_existing_bams() {
    local sample_dirs=(
        "/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_20"
        "/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_21"
        "/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_72"
        "/home/gao/Dropbox/Quote_260203003/Raw_Data/YCL_73"
    )
    
    echo "Checking for existing BAM files..."
    for dir in "${sample_dirs[@]}"; do
        sample_name=$(basename "$dir")
        bam_file="$OUTPUT_DIR/$sample_name.sorted.bam"
        if [ -f "$bam_file" ]; then
            echo "Found BAM for $sample_name: $bam_file"
            analyze_bam_insert_size "$bam_file" "$sample_name"
        else
            echo "No BAM file found for $sample_name"
        fi
    done
}

# Main function
main() {
    echo "Quick Insert Size Analysis Tool"
    echo "================================"
    
    # Check if user provided BAM files as arguments
    if [ $# -gt 0 ]; then
        echo "Analyzing provided BAM files..."
        for bam_file in "$@"; do
            if [ -f "$bam_file" ]; then
                sample_name=$(basename "$bam_file" .bam)
                analyze_bam_insert_size "$bam_file" "$sample_name"
            else
                echo "Warning: BAM file not found: $bam_file"
            fi
        done
    else
        # Check for existing BAM files in results directory
        check_existing_bams
    fi
    
    echo ""
    echo "Analysis complete!"
    echo "Check results in: $OUTPUT_DIR/"
}

# Run main function with all arguments
main "$@"