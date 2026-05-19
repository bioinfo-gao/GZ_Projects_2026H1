#!/usr/bin/env python3
"""
Detailed analysis of primer dimer sequence structure.

This script analyzes merged reads to identify:
1. Read 1 primer regions
2. Read 2 primer regions  
3. Overlapping/complementary regions
4. Random/template regions (if any)
"""

import gzip
import sys
import argparse
from collections import Counter


def reverse_complement(seq):
    """Generate reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))


def find_primer_regions(sequence, max_mismatches=2):
    """Identify primer regions in a sequence."""
    
    # Standard Illumina primer sequences
    primers = {
        'P5_full': 'AATGATACGGCGACCACCGAGATCTACAC',
        'P7_full': 'CAAGCAGAAGACGGCATACGAGAT',
        'Read1_adapter': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
        'Read2_adapter': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
        'Common_adapter': 'AGATCGGAAGAGC'
    }
    
    # Also check reverse complements (since we're looking at merged reads)
    primers_rc = {name + '_RC': reverse_complement(seq) for name, seq in primers.items()}
    all_primers = {**primers, **primers_rc}
    
    found_regions = []
    
    for primer_name, primer_seq in all_primers.items():
        if len(primer_seq) < 8:  # Skip very short sequences
            continue
            
        # Find approximate matches with allowed mismatches
        for i in range(len(sequence) - len(primer_seq) + 1):
            subseq = sequence[i:i+len(primer_seq)]
            mismatches = sum(1 for a, b in zip(subseq, primer_seq) if a != b)
            
            if mismatches <= max_mismatches:
                found_regions.append({
                    'primer': primer_name,
                    'start': i,
                    'end': i + len(primer_seq),
                    'mismatches': mismatches,
                    'matched_seq': subseq,
                    'primer_seq': primer_seq
                })
                break  # Found one match, move to next primer
    
    return found_regions


def analyze_sequence_structure(fastq_file, sample_name, num_reads=100):
    """Analyze the detailed structure of primer dimer sequences."""
    print(f"Analyzing primer dimer structure for {sample_name}...")
    
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
    
    print(f"Analyzing {len(sequences)} sequences...")
    print("")
    
    # Analyze first 10 sequences in detail
    for i, seq in enumerate(sequences[:10]):
        print(f"Sequence {i+1} ({len(seq)} bp):")
        print(f"  Full sequence: {seq}")
        
        # Find primer regions
        primer_regions = find_primer_regions(seq)
        
        if primer_regions:
            print("  Detected primer regions:")
            for region in sorted(primer_regions, key=lambda x: x['start']):
                print(f"    {region['primer']}: positions {region['start']}-{region['end']} "
                      f"({region['mismatches']} mismatches)")
                print(f"      Matched: {region['matched_seq']}")
                print(f"      Primer:  {region['primer_seq']}")
        else:
            print("  No primer regions detected")
        
        # Check for palindromic/self-complementary regions
        # This indicates potential primer-primer interaction
        half_len = len(seq) // 2
        first_half = seq[:half_len]
        second_half_rc = reverse_complement(seq[half_len:])
        
        # Calculate similarity between first half and reverse complement of second half
        min_len = min(len(first_half), len(second_half_rc))
        if min_len > 0:
            matches = sum(1 for a, b in zip(first_half[:min_len], second_half_rc[:min_len]) if a == b)
            similarity = matches / min_len * 100
            print(f"  Self-complementarity: {similarity:.1f}% similarity between "
                  f"first half and RC of second half")
        
        print("")
    
    # Summary statistics
    total_primer_matches = 0
    primer_types = Counter()
    
    for seq in sequences:
        regions = find_primer_regions(seq)
        total_primer_matches += len(regions)
        for region in regions:
            primer_types[region['primer'].split('_')[0]] += 1
    
    print("SUMMARY STATISTICS:")
    print(f"  Total sequences analyzed: {len(sequences)}")
    print(f"  Sequences with primer matches: {sum(1 for seq in sequences if find_primer_regions(seq))}")
    print(f"  Total primer region matches: {total_primer_matches}")
    print("  Primer type distribution:")
    for primer_type, count in primer_types.most_common():
        print(f"    {primer_type}: {count}")
    
    return {
        'total_sequences': len(sequences),
        'sequences_with_primers': sum(1 for seq in sequences if find_primer_regions(seq)),
        'primer_distribution': dict(primer_types)
    }


def main():
    parser = argparse.ArgumentParser(description='Analyze primer dimer sequence structure')
    parser.add_argument('--results-dir', required=True,
                       help='Directory containing fastp results')
    parser.add_argument('--sample', default='YCL_20',
                       help='Sample to analyze (default: YCL_20)')
    
    args = parser.parse_args()
    
    print("PRIMER DIMER SEQUENCE STRUCTURE ANALYSIS")
    print("=" * 60)
    print("")
    print("This analysis identifies:")
    print("- Read 1 primer regions (P5 adapter)")
    print("- Read 2 primer regions (P7 adapter)") 
    print("- Common sequencing adapter regions")
    print("- Self-complementary regions indicating primer-primer interaction")
    print("- Random/template regions (should be minimal in pure primer dimers)")
    print("")
    
    merged_file = f"{args.results_dir}/{args.sample}/{args.sample}_merged.fq.gz"
    results = analyze_sequence_structure(merged_file, args.sample)
    
    if results:
        print("")
        print("INTERPRETATION:")
        print("- Pure primer dimers should show strong P5 and P7 primer signals")
        print("- High self-complementarity indicates primer-primer annealing")
        print("- Minimal random/template regions expected")
        print("- Mixed samples may show both primer and genomic signals")


if __name__ == "__main__":
    main()