#!/usr/bin/env python3
"""
Detailed sequence structure analysis to identify:
- Genomic DNA fragments vs primer sequences
- Adapter positions and boundaries  
- GC-rich regions and their characteristics
- Fragment length patterns
"""

import gzip
import sys
import argparse
from collections import Counter


def reverse_complement(seq):
    """Generate reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))


def calculate_gc_content(sequence):
    """Calculate GC content percentage."""
    if not sequence:
        return 0.0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100


def find_adapter_boundaries(sequence):
    """Find adapter boundaries and classify sequence regions."""
    universal_adapter = 'AGATCGGAAGAGC'
    
    if universal_adapter in sequence:
        adapter_start = sequence.find(universal_adapter)
        genomic_part = sequence[:adapter_start]
        adapter_part = sequence[adapter_start:]
        return {
            'has_adapter': True,
            'genomic_length': len(genomic_part),
            'adapter_length': len(adapter_part),
            'genomic_seq': genomic_part,
            'adapter_seq': adapter_part,
            'boundary_position': adapter_start
        }
    else:
        return {
            'has_adapter': False,
            'genomic_length': len(sequence),
            'adapter_length': 0,
            'genomic_seq': sequence,
            'adapter_seq': '',
            'boundary_position': -1
        }


def analyze_gc_distribution(sequence, window_size=10):
    """Analyze GC content distribution across sequence windows."""
    if len(sequence) < window_size:
        return [calculate_gc_content(sequence)]
    
    gc_windows = []
    for i in range(0, len(sequence) - window_size + 1, window_size//2):
        window = sequence[i:i+window_size]
        gc_windows.append(calculate_gc_content(window))
    
    return gc_windows


def classify_sequence_type(sequence):
    """Classify sequence as genomic-like or primer-like."""
    # Check for repetitive patterns (primer-like)
    repetitive_score = 0
    for i in range(len(sequence) - 3):
        if sequence[i:i+4] in ['AAAA', 'TTTT', 'GGGG', 'CCCC']:
            repetitive_score += 1
    
    # Check complexity
    unique_bases = len(set(sequence))
    complexity_score = unique_bases / len(sequence) if sequence else 0
    
    # Check for known primer patterns
    primer_patterns = ['AGATCGGA', 'CAAGCAGA', 'AATGATAC']
    primer_match = any(pattern in sequence for pattern in primer_patterns)
    
    if repetitive_score > 2 or complexity_score < 0.6 or primer_match:
        return 'primer_like'
    else:
        return 'genomic_like'


def detailed_sequence_analysis(fastq_file, sample_name, num_reads=50):
    """Perform detailed sequence structure analysis."""
    print(f"Detailed sequence structure analysis for {sample_name}...")
    
    sequences = []
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
                    sequences.append(seq)
                    count += 1
    except Exception as e:
        print(f"Error reading {fastq_file}: {e}")
        return None
    
    if not sequences:
        print(f"No sequences found in {fastq_file}")
        return None
    
    print(f"Analyzing {len(sequences)} sequences in detail...\n")
    
    genomic_like_count = 0
    adapter_containing_count = 0
    total_gc_sum = 0
    
    for i, seq in enumerate(sequences):
        print(f"Sequence {i+1} ({len(seq)} bp):")
        print(f"  Full sequence: {seq}")
        
        # Classify sequence type
        seq_type = classify_sequence_type(seq)
        if seq_type == 'genomic_like':
            genomic_like_count += 1
        
        # Find adapter boundaries
        boundary_info = find_adapter_boundaries(seq)
        if boundary_info['has_adapter']:
            adapter_containing_count += 1
            print(f"  Classification: {seq_type.upper()}")
            print(f"  Genomic part: {boundary_info['genomic_seq']} "
                  f"({boundary_info['genomic_length']} bp)")
            print(f"  Adapter part: {boundary_info['adapter_seq']} "
                  f"({boundary_info['adapter_length']} bp)")
        else:
            print(f"  Classification: {seq_type.upper()} (no adapter detected)")
            print(f"  Full sequence is genomic-like fragment")
        
        # GC content analysis
        gc_pct = calculate_gc_content(seq)
        total_gc_sum += gc_pct
        print(f"  GC content: {gc_pct:.1f}%")
        
        # GC distribution analysis
        gc_windows = analyze_gc_distribution(seq)
        if len(gc_windows) > 1:
            gc_variability = max(gc_windows) - min(gc_windows)
            print(f"  GC variability across sequence: {gc_variability:.1f}%")
        
        print("")
    
    # Summary statistics
    mean_gc = total_gc_sum / len(sequences)
    genomic_like_pct = (genomic_like_count / len(sequences)) * 100
    adapter_pct = (adapter_containing_count / len(sequences)) * 100
    
    print("DETAILED ANALYSIS SUMMARY:")
    print(f"  Total sequences analyzed: {len(sequences)}")
    print(f"  Genomic-like sequences: {genomic_like_count} ({genomic_like_pct:.1f}%)")
    print(f"  Adapter-containing sequences: {adapter_containing_count} ({adapter_pct:.1f}%)")
    print(f"  Mean GC content: {mean_gc:.1f}%")
    
    if genomic_like_pct > 80:
        print("\n✅ CONCLUSION: Sequences are predominantly GENOMIC-LIKE")
        print("   This strongly supports DNA degradation rather than primer dimers.")
    else:
        print("\n⚠️  CONCLUSION: Significant primer-like sequences detected")
        print("   Consider mixed contamination or technical artifacts.")
    
    return {
        'genomic_like_pct': genomic_like_pct,
        'adapter_pct': adapter_pct,
        'mean_gc': mean_gc
    }


def main():
    parser = argparse.ArgumentParser(description='Detailed sequence structure analysis')
    parser.add_argument('--results-dir', required=True,
                       help='Directory containing fastp results')
    parser.add_argument('--sample', default='YCL_20',
                       help='Sample to analyze in detail (default: YCL_20)')
    
    args = parser.parse_args()
    
    print("DETAILED SEQUENCE STRUCTURE ANALYSIS")
    print("=" * 50)
    print("")
    print("This analysis examines individual sequences to determine:")
    print("- Whether they are genomic-like or primer-like")
    print("- Adapter boundary positions")  
    print("- GC content distribution patterns")
    print("- Sequence complexity and repetitiveness")
    print("")
    
    merged_file = f"{args.results_dir}/{args.sample}/{args.sample}_merged.fq.gz"
    results = detailed_sequence_analysis(merged_file, args.sample)
    
    if results:
        print("\nFor comprehensive analysis across all samples, run:")
        print("python3 analyze_degradation_evidence.py --results-dir results_fastp")


if __name__ == "__main__":
    main()