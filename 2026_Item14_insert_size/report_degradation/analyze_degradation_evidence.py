#!/usr/bin/env python3
"""
Comprehensive DNA degradation analysis to distinguish between 
primer dimers and degraded genomic DNA.

This script analyzes:
1. Sequence composition and structure
2. GC content bias patterns  
3. Adapter contamination levels
4. Fragment length distributions
5. Sequence complexity metrics
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


def find_adapter_regions(sequence):
    """Find Illumina adapter regions in sequence."""
    adapters = {
        'universal_adapter': 'AGATCGGAAGAGC',
        'read1_adapter': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
        'read2_adapter': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    }
    
    found_adapters = []
    for adapter_name, adapter_seq in adapters.items():
        if adapter_seq in sequence:
            start_pos = sequence.find(adapter_seq)
            found_adapters.append({
                'name': adapter_name,
                'start': start_pos,
                'end': start_pos + len(adapter_seq),
                'sequence': adapter_seq
            })
    
    return found_adapters


def analyze_sequence_structure(fastq_file, sample_name, num_reads=1000):
    """Analyze sequence structure for degradation evidence."""
    print(f"Analyzing degradation evidence for {sample_name}...")
    
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
    
    total_sequences = len(sequences)
    
    # 1. GC content analysis
    gc_contents = [calculate_gc_content(seq) for seq in sequences]
    mean_gc = sum(gc_contents) / len(gc_contents)
    high_gc_count = sum(1 for gc in gc_contents if gc > 60)
    
    # 2. Adapter contamination analysis
    adapter_contamination = 0
    adapter_details = []
    for seq in sequences[:100]:  # Analyze first 100 for detailed adapter analysis
        adapters = find_adapter_regions(seq)
        if adapters:
            adapter_contamination += 1
            adapter_details.extend(adapters)
    
    adapter_pct = (adapter_contamination / 100) * 100 if sequences else 0
    
    # 3. Sequence uniqueness
    unique_sequences = len(set(sequences))
    uniqueness_pct = (unique_sequences / total_sequences) * 100
    
    # 4. Length distribution
    lengths = [len(seq) for seq in sequences]
    mean_length = sum(lengths) / len(lengths)
    length_std = (sum((l - mean_length) ** 2 for l in lengths) / len(lengths)) ** 0.5
    
    # 5. Complexity analysis (k-mer diversity)
    kmer_counter = Counter()
    for seq in sequences[:500]:  # Analyze first 500 sequences
        for i in range(len(seq) - 3):
            kmer = seq[i:i+4]
            kmer_counter[kmer] += 1
    
    total_kmers = sum(kmer_counter.values())
    unique_kmers = len(kmer_counter)
    kmer_diversity = unique_kmers / total_kmers if total_kmers > 0 else 0
    
    # Print results
    print(f"  Total sequences analyzed: {total_sequences}")
    print(f"  Mean GC content: {mean_gc:.1f}%")
    print(f"  Sequences with GC > 60%: {high_gc_count} ({high_gc_count/total_sequences*100:.1f}%)")
    print(f"  Adapter contamination: {adapter_pct:.1f}% ({adapter_contamination}/100)")
    print(f"  Sequence uniqueness: {uniqueness_pct:.1f}% ({unique_sequences} unique)")
    print(f"  Mean fragment length: {mean_length:.1f} bp (SD: {length_std:.1f})")
    print(f"  4-mer diversity: {kmer_diversity:.4f}")
    
    # Show example sequences
    print(f"\n  Example sequences:")
    for i, seq in enumerate(sequences[:5]):
        gc_pct = calculate_gc_content(seq)
        adapters = find_adapter_regions(seq)
        adapter_info = f" + adapter" if adapters else ""
        print(f"    {i+1}. {seq[:50]}{'...' if len(seq) > 50 else ''} "
              f"({len(seq)} bp, {gc_pct:.1f}% GC{adapter_info})")
    
    return {
        'mean_gc': mean_gc,
        'adapter_pct': adapter_pct,
        'uniqueness_pct': uniqueness_pct,
        'mean_length': mean_length,
        'kmer_diversity': kmer_diversity,
        'high_gc_pct': high_gc_count/total_sequences*100
    }


def main():
    parser = argparse.ArgumentParser(description='Analyze DNA degradation evidence')
    parser.add_argument('--results-dir', required=True,
                       help='Directory containing fastp results')
    
    args = parser.parse_args()
    
    samples = ['YCL_20', 'YCL_21', 'YCL_72', 'YCL_73']
    
    print("DNA DEGRADATION EVIDENCE ANALYSIS")
    print("=" * 60)
    print("")
    print("This analysis distinguishes between:")
    print("- Primer dimers (technical artifacts)")
    print("- Degraded genomic DNA (biological reality)")
    print("")
    print("Key discriminators:")
    print("- Sequence uniqueness (>95% suggests genomic DNA)")
    print("- Adapter contamination (<10% suggests real fragments)")  
    print("- GC bias pattern (moderate 60-65% vs extreme >70%)")
    print("- Sequence structure (genomic-like vs primer-like)")
    print("")
    
    all_results = {}
    for sample in samples:
        merged_file = f"{args.results_dir}/{sample}/{sample}_merged.fq.gz"
        results = analyze_sequence_structure(merged_file, sample)
        if results:
            all_results[sample] = results
        print("")
    
    # Summary table
    if all_results:
        print("DEGRADATION EVIDENCE SUMMARY")
        print("-" * 80)
        print(f"{'Sample':<15} {'GC%':<8} {'Adapter%':<12} {'Unique%':<12} {'Length':<10} {'k-mer Div'}")
        print("-" * 80)
        for sample, results in all_results.items():
            print(f"{sample:<15} {results['mean_gc']:<8.1f} {results['adapter_pct']:<12.1f} "
                  f"{results['uniqueness_pct']:<12.1f} {results['mean_length']:<10.1f} {results['kmer_diversity']:.4f}")
        
        print("")
        print("INTERPRETATION:")
        print("- High uniqueness (>95%) + low adapter (<10%) = Degraded genomic DNA")
        print("- Low uniqueness (<50%) + high adapter (>50%) = Primer dimers")
        print("- Moderate GC (60-65%) = Selective AT degradation")
        print("- Extreme GC (>70%) = Primer design bias")
        print("")
        
        # Classification
        degradation_samples = []
        for sample, results in all_results.items():
            if (results['uniqueness_pct'] > 95 and 
                results['adapter_pct'] < 10 and 
                60 <= results['mean_gc'] <= 70):
                degradation_samples.append(sample)
        
        if len(degradation_samples) == len(samples):
            print("✅ CONCLUSION: ALL SAMPLES SHOW EVIDENCE OF DNA DEGRADATION")
            print("   Not primer dimers - these are real but degraded DNA fragments.")
            print("   Root cause likely sample quality issues before library prep.")
        elif len(degradation_samples) > 0:
            print(f"✅ CONCLUSION: {len(degradation_samples)} samples show DNA degradation")
            print("   Remaining samples may have mixed contamination.")
        else:
            print("⚠️  CONCLUSION: Evidence suggests primer dimer contamination")
            print("   Consider technical artifact rather than biological degradation.")


if __name__ == "__main__":
    main()