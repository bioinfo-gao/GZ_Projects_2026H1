#!/usr/bin/env python3
"""
GC Content Analysis for Primer Dimer Detection

This script analyzes GC content distribution to confirm primer dimer contamination.
Mammalian genomic DNA should have GC content of 40-45%, with maximum rarely exceeding 55%.
Observed GC content of 60-65% strongly indicates technical artifact.
"""

import gzip
import sys
import argparse
from collections import Counter


def calculate_gc_content(sequence):
    """Calculate GC content percentage for a sequence."""
    if not sequence:
        return 0.0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100


def analyze_gc_distribution(fastq_file, sample_name, num_reads=1000):
    """Analyze GC content distribution in merged reads."""
    print(f"Analyzing GC content for {sample_name}...")
    
    gc_contents = []
    lengths = []
    
    try:
        with gzip.open(fastq_file, 'rt') as f:
            count = 0
            while count < num_reads:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip().upper()
                f.readline()  # plus line
                f.readline()  # quality line
                
                if seq:
                    gc_pct = calculate_gc_content(seq)
                    gc_contents.append(gc_pct)
                    lengths.append(len(seq))
                    count += 1
    except Exception as e:
        print(f"Error reading {fastq_file}: {e}")
        return None
    
    if not gc_contents:
        print(f"No sequences found in {fastq_file}")
        return None
    
    # Calculate statistics
    mean_gc = sum(gc_contents) / len(gc_contents)
    median_gc = sorted(gc_contents)[len(gc_contents)//2]
    min_gc = min(gc_contents)
    max_gc = max(gc_contents)
    
    # Count reads by GC ranges
    normal_gc = sum(1 for gc in gc_contents if gc <= 55)
    high_gc = sum(1 for gc in gc_contents if gc > 55)
    very_high_gc = sum(1 for gc in gc_contents if gc > 60)
    
    print(f"  Total reads analyzed: {len(gc_contents)}")
    print(f"  Mean GC content: {mean_gc:.1f}%")
    print(f"  Median GC content: {median_gc:.1f}%")
    print(f"  GC range: {min_gc:.1f}% - {max_gc:.1f}%")
    print(f"  Reads with GC <= 55%: {normal_gc} ({normal_gc/len(gc_contents)*100:.1f}%)")
    print(f"  Reads with GC > 55%: {high_gc} ({high_gc/len(gc_contents)*100:.1f}%)")
    print(f"  Reads with GC > 60%: {very_high_gc} ({very_high_gc/len(gc_contents)*100:.1f}%)")
    
    return {
        'mean_gc': mean_gc,
        'median_gc': median_gc,
        'high_gc_pct': high_gc/len(gc_contents)*100,
        'very_high_gc_pct': very_high_gc/len(gc_contents)*100
    }


def main():
    parser = argparse.ArgumentParser(description='Analyze GC content for primer dimer detection')
    parser.add_argument('--results-dir', required=True,
                       help='Directory containing fastp results')
    
    args = parser.parse_args()
    
    samples = ['YCL_20', 'YCL_21', 'YCL_72', 'YCL_73']
    
    print("GC CONTENT ANALYSIS FOR PRIMER DIMER DETECTION")
    print("=" * 60)
    print("")
    print("Background: Normal mammalian genomic DNA has GC content of 40-45%")
    print("Maximum observed in real samples: typically < 55%")
    print("Primer dimers often show elevated GC (>60%) due to primer design")
    print("")
    
    all_results = {}
    for sample in samples:
        merged_file = f"{args.results-dir}/{sample}/{sample}_merged.fq.gz"
        results = analyze_gc_distribution(merged_file, sample)
        if results:
            all_results[sample] = results
        print("")
    
    # Summary table
    if all_results:
        print("SUMMARY TABLE")
        print("-" * 70)
        print(f"{'Sample':<15} {'Mean GC%':<12} {'>55%':<12} {'>60%'}")
        print("-" * 70)
        for sample, results in all_results.items():
            mean_gc = results['mean_gc']
            high_gc = results['high_gc_pct']
            very_high_gc = results['very_high_gc_pct']
            print(f"{sample:<15} {mean_gc:<12.1f} {high_gc:<12.1f} {very_high_gc:.1f}")
        
        print("")
        print("INTERPRETATION:")
        print("- Normal mammalian DNA: GC% ≤ 55%")
        print("- Primer dimers/artifacts: Often GC% > 60%")
        print("- Consistent high GC across samples suggests technical artifact")
        print("")
        print("CONCLUSION:")
        if all(results['mean_gc'] > 55 for results in all_results.values()):
            print("⚠️  ALL SAMPLES SHOW ABNORMALLY HIGH GC CONTENT!")
            print("   This strongly supports primer dimer contamination.")
            print("   Data is likely unsuitable for intended biological analysis.")
        else:
            print("GC content appears within normal biological ranges.")


if __name__ == "__main__":
    main()