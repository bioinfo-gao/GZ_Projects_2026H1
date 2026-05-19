#!/usr/bin/env python3
"""
Detailed sequence complexity analysis to distinguish between 
genomic DNA and primer dimers.

This script analyzes:
1. Sequence uniqueness (duplicate rate)
2. k-mer diversity 
3. Base composition bias
4. Sequence clustering patterns
"""

import gzip
import sys
import argparse
from collections import Counter, defaultdict


def analyze_sequence_complexity(fastq_file, sample_name, num_reads=5000):
    """Analyze sequence complexity in detail."""
    print(f"Analyzing sequence complexity for {sample_name}...")
    
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
    
    total_reads = len(sequences)
    
    # 1. Duplicate rate analysis
    seq_counter = Counter(sequences)
    unique_sequences = len(seq_counter)
    duplicate_rate = (total_reads - unique_sequences) / total_reads * 100
    
    # Top duplicates
    top_duplicates = seq_counter.most_common(10)
    
    # 2. Analyze first 50bp (as suggested by user)
    first_50bp_sequences = [seq[:50] for seq in sequences if len(seq) >= 50]
    if first_50bp_sequences:
        first_50bp_counter = Counter(first_50bp_sequences)
        unique_50bp = len(first_50bp_counter)
        duplicate_rate_50bp = (len(first_50bp_sequences) - unique_50bp) / len(first_50bp_sequences) * 100
        top_50bp_duplicates = first_50bp_counter.most_common(10)
    else:
        duplicate_rate_50bp = 0
        top_50bp_duplicates = []
    
    # 3. k-mer analysis (k=4)
    kmer_counter = Counter()
    for seq in sequences[:1000]:  # Analyze first 1000 sequences
        for i in range(len(seq) - 3):
            kmer = seq[i:i+4]
            kmer_counter[kmer] += 1
    
    total_kmers = sum(kmer_counter.values())
    unique_kmers = len(kmer_counter)
    kmer_diversity = unique_kmers / total_kmers if total_kmers > 0 else 0
    
    # Most frequent 4-mers
    top_kmers = kmer_counter.most_common(10)
    
    # 4. Base composition analysis
    base_counter = Counter()
    for seq in sequences[:1000]:
        base_counter.update(seq)
    
    total_bases = sum(base_counter.values())
    base_freq = {base: count/total_bases*100 for base, count in base_counter.items()}
    
    print(f"  Total reads analyzed: {total_reads}")
    print(f"  Unique sequences: {unique_sequences} ({unique_sequences/total_reads*100:.1f}%)")
    print(f"  Duplicate rate: {duplicate_rate:.1f}%")
    print(f"  First 50bp unique: {unique_50bp if first_50bp_sequences else 'N/A'}")
    if first_50bp_sequences:
        print(f"  First 50bp duplicate rate: {duplicate_rate_50bp:.1f}%")
    print(f"  k-mer diversity (k=4): {kmer_diversity:.4f}")
    print(f"  Base frequencies: A={base_freq.get('A',0):.1f}%, C={base_freq.get('C',0):.1f}%, G={base_freq.get('G',0):.1f}%, T={base_freq.get('T',0):.1f}%")
    
    # Show top duplicates if they exist
    if top_duplicates and top_duplicates[0][1] > 1:
        print(f"  Top duplicate sequences:")
        for i, (seq, count) in enumerate(top_duplicates[:5]):
            if count > 1:
                print(f"    {i+1}. {seq[:30]}... ({count} times)")
    
    if top_50bp_duplicates and top_50bp_duplicates[0][1] > 1:
        print(f"  Top 50bp duplicate sequences:")
        for i, (seq, count) in enumerate(top_50bp_duplicates[:5]):
            if count > 1:
                print(f"    {i+1}. {seq[:30]}... ({count} times)")
    
    return {
        'duplicate_rate': duplicate_rate,
        'duplicate_rate_50bp': duplicate_rate_50bp if first_50bp_sequences else 0,
        'kmer_diversity': kmer_diversity,
        'unique_sequences_pct': unique_sequences/total_reads*100,
        'top_duplicates': top_duplicates[:5] if top_duplicates else []
    }


def main():
    parser = argparse.ArgumentParser(description='Analyze sequence complexity for primer dimer detection')
    parser.add_argument('--results-dir', required=True,
                       help='Directory containing fastp results')
    
    args = parser.parse_args()
    
    samples = ['YCL_20', 'YCL_21', 'YCL_72', 'YCL_73']
    
    print("DETAILED SEQUENCE COMPLEXITY ANALYSIS")
    print("=" * 60)
    print("")
    print("Analysis focuses on:")
    print("- Sequence duplication rates (primer dimers show high duplication)")
    print("- First 50bp sequence analysis (as suggested by user)")
    print("- k-mer diversity (low diversity suggests technical artifacts)")
    print("- Base composition bias")
    print("")
    
    all_results = {}
    for sample in samples:
        merged_file = f"{args.results_dir}/{sample}/{sample}_merged.fq.gz"
        results = analyze_sequence_complexity(merged_file, sample)
        if results:
            all_results[sample] = results
        print("")
    
    # Summary table
    if all_results:
        print("COMPLEXITY SUMMARY")
        print("-" * 80)
        print(f"{'Sample':<15} {'Dup%':<10} {'50bp Dup%':<12} {'k-mer Div':<12} {'Unique%'}")
        print("-" * 80)
        for sample, results in all_results.items():
            dup_rate = results['duplicate_rate']
            dup_50bp = results['duplicate_rate_50bp']
            kmer_div = results['kmer_diversity']
            unique_pct = results['unique_sequences_pct']
            print(f"{sample:<15} {dup_rate:<10.1f} {dup_50bp:<12.1f} {kmer_div:<12.4f} {unique_pct:.1f}")
        
        print("")
        print("INTERPRETATION:")
        print("- Genomic DNA: Low duplication (<10%), high k-mer diversity (>0.1)")
        print("- Primer dimers: High duplication (>50%), low k-mer diversity (<0.05)")
        print("- First 50bp analysis is critical for short fragments")
        print("")
        
        # Check for primer dimer signatures
        primer_dimer_samples = []
        for sample, results in all_results.items():
            if (results['duplicate_rate'] > 30 or 
                results['duplicate_rate_50bp'] > 30 or 
                results['kmer_diversity'] < 0.05):
                primer_dimer_samples.append(sample)
        
        if primer_dimer_samples:
            print("⚠️  PRIMER DIMER SIGNATURES DETECTED IN:")
            for sample in primer_dimer_samples:
                print(f"   - {sample}")
        else:
            print("✅ No strong primer dimer signatures detected")
            print("   However, high GC content (60-65%) remains concerning")


if __name__ == "__main__":
    main()