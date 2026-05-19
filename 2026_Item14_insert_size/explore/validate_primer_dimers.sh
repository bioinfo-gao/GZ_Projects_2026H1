#!/bin/bash

# Validate if the data contains primer dimers
# This script checks for adapter/primer sequences in merged reads

set -e

echo "Validating Primer Dimer Contamination..."
echo "======================================="

# Common Illumina adapter sequences
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

OUTPUT_DIR="/home/gao/projects/2026_Item14_insert_size/validation_results"
mkdir -p "$OUTPUT_DIR"

SAMPLES=("YCL_20" "YCL_21" "YCL_72" "YCL_73")

for sample in "${SAMPLES[@]}"; do
    echo "Analyzing $sample..."
    
    # Get merged reads file
    MERGED_FILE="/home/gao/projects/2026_Item14_insert_size/results_fastp/$sample/${sample}_merged.fq.gz"
    
    if [ ! -f "$MERGED_FILE" ]; then
        echo "  Warning: Merged file not found for $sample"
        continue
    fi
    
    # Extract first 1000 reads for analysis
    TEMP_FASTA="$OUTPUT_DIR/${sample}_test_reads.fa"
    zcat "$MERGED_FILE" | head -4000 | awk 'NR%4==1 {print ">" substr($0,2)} NR%4==2 {print}' > "$TEMP_FASTA"
    
    # Check for adapter sequences at beginning of reads
    echo "  Checking for adapter contamination..."
    ADAPTER_COUNT_R1=$(grep -c "^$ADAPTER_R1" "$TEMP_FASTA")
    ADAPTER_COUNT_R2=$(grep -c "^$ADAPTER_R2" "$TEMP_FASTA")
    
    echo "  Adapter R1 matches: $ADAPTER_COUNT_R1"
    echo "  Adapter R2 matches: $ADAPTER_COUNT_R2"
    
    # Check read length distribution
    echo "  Read length distribution (first 1000 reads):"
    awk 'NR%2==0 {print length}' "$TEMP_FASTA" | sort -n | uniq -c > "$OUTPUT_DIR/${sample}_length_dist.txt"
    head -5 "$OUTPUT_DIR/${sample}_length_dist.txt"
    
    # Calculate percentage of reads with adapter
    TOTAL_READS=1000
    ADAPTER_TOTAL=$((ADAPTER_COUNT_R1 + ADAPTER_COUNT_R2))
    ADAPTER_PERCENTAGE=$(echo "scale=2; $ADAPTER_TOTAL * 100 / $TOTAL_READS" | bc)
    echo "  Adapter contamination: ${ADAPTER_PERCENTAGE}%"
    
    echo ""
done

echo "Validation complete!"
echo "Check detailed results in: $OUTPUT_DIR/"
echo ""
echo "Interpretation guide:"
echo "- If adapter contamination > 50%: Likely primer dimers"
echo "- If most reads are 50-70 bp: Consistent with primer dimers"  
echo "- If adapter sequences dominate: Confirms technical artifact"