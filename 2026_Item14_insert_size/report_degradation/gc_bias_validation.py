#!/usr/bin/env python3
"""
GC bias validation for DNA degradation analysis.

This script validates that the observed GC bias (60-65%) is consistent 
with DNA degradation rather than technical artifacts like primer dimers.
"""

import gzip
import sys
import argparse
import matplotlib.pyplot as plt
from collections import Counter


def calculate_gc_content(sequence):
    """Calculate GC content percentage."""
    if not sequence:
        return 0.0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100


def analyze_gc_distribution(fastq_file, sample_name, num_reads=5000):
    """Analyze GC content distribution across sequences."""
    print(f"Analyzing GC distribution for {sample_name}...")
    
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
        return None, None
    
    if not gc_contents:
        print(f"No sequences found in {fastq_file}")
        return None, None
    
    return gc_contents, lengths


def validate_gc_bias_pattern(gc_contents):
    """Validate GC bias pattern for degradation vs primer dimers."""
    
    mean_gc = sum(gc_contents) / len(gc_contents)
    median_gc = sorted(gc_contents)[len(gc_contents)//2]
    
    # Count sequences in different GC ranges
    gc_ranges = {
        'low': sum(1 for gc in gc_contents if gc <= 45),      # Normal genomic range
        'moderate': sum(1 for gc in gc_contents if 46 <= gc <= 55),  # Slightly elevated
        'high': sum(1 for gc in gc_contents if 56 <= gc <= 70),     # Degradation range  
        'very_high': sum(1 for gc in gc_contents if gc > 70)        # Primer dimer range
    }
    
    total = len(gc_contents)
    percentages = {k: v/total*100 for k, v in gc_ranges.items()}
    
    print(f"  Mean GC: {mean_gc:.1f}%")
    print(f"  Median GC: {median_gc:.1f}%")
    print(f"  GC Distribution:")
    print(f"    ≤45% (normal): {gc_ranges['low']} ({percentages['low']:.1f}%)")
    print(f"    46-55% (elevated): {gc_ranges['moderate']} ({percentages['moderate']:.1f}%)")
    print(f"    56-70% (degraded): {gc_ranges['high']} ({percentages['high']:.1f}%)")
    print(f"    >70% (primer-like): {gc_ranges['very_high']} ({percentages['very_high']:.1f}%)")
    
    # Classification logic
    if percentages['high'] > 40 and percentages['very_high'] < 20:
        classification = "DNA DEGRADATION"
        evidence = "Moderate GC enrichment (56-70%) consistent with AT-rich region loss"
    elif percentages['very_high'] > 30:
        classification = "PRIMER DIMER CONTAMINATION"
        evidence = "Extreme GC enrichment (>70%) consistent with primer design bias"
    else:
        classification = "MIXED OR UNCERTAIN"
        evidence = "Insufficient GC bias pattern for clear classification"
    
    print(f"  Classification: {classification}")
    print(f"  Evidence: {evidence}")
    
    return {
        'mean_gc': mean_gc,
        'median_gc': median_gc,
        'percentages': percentages,
        'classification': classification,
        'evidence': evidence
    }


def main():
    parser = argparse.ArgumentParser(description='Validate GC bias for DNA degradation')
    parser.add_argument('--results-dir', required=True,
                       help='Directory containing fastp results')
    
    args = parser.parse_args()
    
    samples = ['YCL_20', 'YCL_21', 'YCL_72', 'YCL_73']
    
    print("GC BIAS VALIDATION FOR DNA DEGRADATION ANALYSIS")
    print("=" * 60)
    print("")
    print("This analysis validates whether high GC content is consistent with:")
    print("- DNA degradation (moderate GC enrichment: 56-70%)")
    print("- Primer dimer contamination (extreme GC enrichment: >70%)")
    print("")
    
    all_results = {}
    for sample in samples:
        merged_file = f"{args.results_dir}/{sample}/{sample}_merged.fq.gz"
        gc_contents, lengths = analyze_gc_distribution(merged_file, sample)
        if gc_contents:
            results = validate_gc_bias_pattern(gc_contents)
            all_results[sample] = results
        print("")
    
    # Summary table
    if all_results:
        print("GC BIAS VALIDATION SUMMARY")
        print("-" * 80)
        print(f"{'Sample':<15} {'Mean GC%':<12} {'56-70%':<12} {'>70%':<12} {'Classification'}")
        print("-" * 80)
        for sample, results in all_results.items():
            mean_gc = results['mean_gc']
            degraded_pct = results['percentages']['high']
            primer_pct = results['percentages']['very_high']
            classification = results['classification']
            print(f"{sample:<15} {mean_gc:<12.1f} {degraded_pct:<12.1f} {primer_pct:<12.1f} {classification}")
        
        print("")
        degradation_samples = [s for s, r in all_results.items() if r['classification'] == 'DNA DEGRADATION']
        primer_samples = [s for s, r in all_results.items() if r['classification'] == 'PRIMER DIMER CONTAMINATION']
        
        if len(degradation_samples) == len(samples):
            print("✅ CONCLUSION: ALL SAMPLES SHOW GC BIAS CONSISTENT WITH DNA DEGRADATION")
            print("   Moderate GC enrichment (56-70%) supports selective AT-rich region loss.")
            print("   This is biological degradation, not technical primer dimer artifact.")
        elif len(primer_samples) == len(samples):
            print("⚠️  CONCLUSION: ALL SAMPLES SHOW GC BIAS CONSISTENT WITH PRIMER DIMERS")
            print("   Extreme GC enrichment (>70%) indicates technical contamination.")
        else:
            print("❓ CONCLUSION: MIXED PATTERNS DETECTED")
            if degradation_samples:
                print(f"   Degradation samples: {', '.join(degradation_samples)}")
            if primer_samples:
                print(f"   Primer dimer samples: {', '.join(primer_samples)}")


if __name__ == "__main__":
    main()