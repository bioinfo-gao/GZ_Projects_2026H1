#!/usr/bin/env python3
"""
Analyze the first 50bp of merged reads to determine if they represent 
primer dimers or genuine genomic sequences.

This script extracts the first 50bp (or full length if shorter) from 
merged reads and analyzes:
1. Sequence uniqueness/repetition rates
2. k-mer frequency distribution  
3. Sequence complexity metrics
4. Common sequence patterns
"""

import gzip
import sys
import argparse
from collections import Counter


def analyze_sequence_complexity(sequences, sample_name):
    """Analyze sequence complexity and repetition patterns."""
    print(f"Analyzing first 50bp sequences for {sample_name}...")
    
    total_sequences = len(sequences)
    if total_sequences == 0:
        print("  No sequences to analyze")
        return
    
    # 1. Uniqueness analysis
    unique_sequences = set(sequences)
    unique_count = len(unique_sequences)
    uniqueness_rate = unique_count / total_sequences * 100
    
    # Top repeated sequences
    seq_counter = Counter(sequences)
    top_repeated = seq_counter.most_common(10)
    
    print(f"  Total sequences analyzed: {total_sequences}")
    print(f"  Unique sequences: {unique_count} ({uniqueness_rate:.1f}%)")
    print(f"  Most repeated sequences:")
    for i, (seq, count) in enumerate(top_repeated[:5], 1):
        percentage = count / total_sequences * 100
        print(f"    {i}. '{seq[:30]}{'...' if len(seq) > 30 else ''}' - {count} times ({percentage:.1f}%)")
    
    # 2. Length distribution of extracted sequences
    lengths = [len(seq) for seq in sequences]
    mean_length = sum(lengths) / len(lengths)
    print(f"  Mean extracted length: {mean_length:.1f} bp")
    
    # 3. Complexity analysis using k-mers (k=4)
    all_4mers = []
    for seq in sequences:
        if len(seq) >= 4:
            for i in range(len(seq) - 3):
                all_4mers.append(seq[i:i+4])
    
    if all_4mers:
        unique_4mers = len(set(all_4mers))
        total_4mers = len(all_4mers)
        kmer_diversity = unique_4mers / total_4mers * 100
        print(f"  4-mer diversity: {unique_4mers}/{total_4mers} ({kmer_diversity:.1f}%)")
    else:
        kmer_diversity = 0
        print("  Insufficient sequence length for k-mer analysis")
    
    # 4. Check for common primer/adapter patterns
    common_patterns = {
        'Illumina_R1': 'AGATCGGAAGAGC',
        'Illumina_R2': 'AGATCGGAAGAGCGTC',
        'Poly_A': 'AAAAAAAA',
        'Poly_T': 'TTTTTTTT',
        'Poly_G': 'GGGGGGGG',  
        'Poly_C': 'CCCCCCCC'
    }
    
    pattern_matches = {}
    for pattern_name, pattern_seq in common_patterns.items():
        matches = sum(1 for seq in sequences if pattern_seq in seq)
        if matches > 0:
            pattern_matches[pattern_name] = matches
    
    if pattern_matches:
        print("  Detected common patterns:")
        for pattern, count in pattern_matches.items():
            pct = count / total_sequences * 100
            print(f"    {pattern}: {count} ({pct:.1f}%)")
    
    return {
        'uniqueness_rate': uniqueness_rate,
        'kmer_diversity': kmer_diversity,
        'top_repeated': top_repeated[:3],
        'pattern_matches': pattern_matches
    }


def extract_first_50bp(fastq_file, max_reads=5000):
    """Extract first 50bp from merged reads."""
    sequences = []
    try:
        with gzip.open(fastq_file, 'rt') as f:
            count = 0
            while count < max_reads:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip().upper()
                f.readline()  # plus line
                f.readline()  # quality line
                
                if seq:
                    # Extract first 50bp or full sequence if shorter
                    extracted = seq[:50]
                    sequences.append(extracted)
                    count += 1
    except Exception as e:
        print(f"Error reading {fastq_file}: {e}")
        return None
    
    return sequences


def main():
    parser = argparse.ArgumentParser(description='Analyze first 50bp of merged reads for primer dimer detection')
    parser.add_argument('--results-dir', required=True,
                       help='Directory containing fastp results')
    parser.add_argument('--max-reads', type=int, default=5000,
                       help='Maximum number of reads to analyze per sample')
    
    args = parser.parse_args()
    
    samples = ['YCL_20', 'YCL_21', 'YCL_72', 'YCL_73']
    
    print("FIRST 50BP SEQUENCE COMPLEXITY ANALYSIS")
    print("=" * 60)
    print("")
    print("Primer dimer characteristics:")
    print("- Low sequence uniqueness (<10% unique sequences)")
    print("- High repetition of identical sequences")
    print("- Low k-mer diversity (<50%)")
    print("- Presence of adapter/primer patterns")
    print("")
    
    all_results = {}
    for sample in samples:
        merged_file = f"{args.results_dir}/{sample}/{sample}_merged.fq.gz"
        sequences = extract_first_50bp(merged_file, args.max_reads)
        if sequences:
            results = analyze_sequence_complexity(sequences, sample)
            all_results[sample] = results
        else:
            print(f"Failed to extract sequences for {sample}")
        print("")
    
    # Summary interpretation
    if all_results:
        print("INTERPRETATION SUMMARY")
        print("-" * 40)
        
        primer_dimer_evidence = []
        for sample, results in all_results.items():
            evidence = []
            
            # Check uniqueness
            if results['uniqueness_rate'] < 10:
                evidence.append("low_uniqueness")
            
            # Check k-mer diversity  
            if results['kmer_diversity'] < 50:
                evidence.append("low_kmer_diversity")
                
            # Check pattern matches
            if results['pattern_matches']:
                evidence.append("adapter_patterns")
            
            # Check top repetition
            if results['top_repeated']:
                top_count = results['top_repeated'][0][1]
                if top_count > len(sequences) * 0.05:  # Top sequence >5% of total
                    evidence.append("high_repetition")
            
            if evidence:
                primer_dimer_evidence.append(sample)
                print(f"{sample}: Evidence for primer dimers - {', '.join(evidence)}")
            else:
                print(f"{sample}: No strong evidence for primer dimers")
        
        print("")
        if len(primer_dimer_evidence) == len(samples):
            print("⚠️  CONCLUSION: ALL SAMPLES SHOW STRONG EVIDENCE OF PRIMER DIMER CONTAMINATION")
            print("   This is consistent with the abnormally high GC content (60-65%)")
            print("   and small insert sizes (58-60 bp) observed previously.")
        elif len(primer_dimer_evidence) > 0:
            print(f"⚠️  CONCLUSION: {len(primer_dimer_evidence)} samples show primer dimer evidence")
        else:
            print("✅ CONCLUSION: No strong evidence of primer dimer contamination detected")
            print("   However, abnormally high GC content (60-65%) remains concerning")


if __name__ == "__main__":
    main()