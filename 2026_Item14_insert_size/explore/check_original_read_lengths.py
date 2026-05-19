#!/usr/bin/env python3
"""
Check original FASTQ read lengths to verify if they are short (40-60 bp)
as expected for primer dimers or small fragments.
"""

import gzip
import sys


def check_fastq_lengths(fastq_file, max_reads=1000):
    """Check read lengths in a FASTQ file."""
    lengths = []
    try:
        with gzip.open(fastq_file, 'rt') as f:
            count = 0
            while count < max_reads:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()  # plus line
                f.readline()  # quality line
                
                if seq:
                    lengths.append(len(seq))
                    count += 1
    except Exception as e:
        print(f"Error reading {fastq_file}: {e}")
        return None
    
    if not lengths:
        print(f"No sequences found in {fastq_file}")
        return None
    
    # Calculate statistics
    min_len = min(lengths)
    max_len = max(lengths)
    mean_len = sum(lengths) / len(lengths)
    
    # Count reads by length ranges
    short_reads = sum(1 for l in lengths if l <= 60)
    medium_reads = sum(1 for l in lengths if 61 <= l <= 100)
    long_reads = sum(1 for l in lengths if l > 100)
    
    print(f"Read length analysis for {fastq_file}:")
    print(f"  Total reads analyzed: {len(lengths)}")
    print(f"  Length range: {min_len}-{max_len} bp")
    print(f"  Mean length: {mean_len:.1f} bp")
    print(f"  Reads <= 60 bp: {short_reads} ({short_reads/len(lengths)*100:.1f}%)")
    print(f"  Reads 61-100 bp: {medium_reads} ({medium_reads/len(lengths)*100:.1f}%)")
    print(f"  Reads > 100 bp: {long_reads} ({long_reads/len(lengths)*100:.1f}%)")
    
    return {
        'min': min_len,
        'max': max_len,
        'mean': mean_len,
        'short_pct': short_reads/len(lengths)*100,
        'total': len(lengths)
    }


def main():
    samples = ['YCL_20', 'YCL_21', 'YCL_72', 'YCL_73']
    base_dir = '/home/gao/Dropbox/Quote_260203003/Raw_Data'
    
    print("ORIGINAL FASTQ READ LENGTH ANALYSIS")
    print("=" * 50)
    print("")
    
    all_results = {}
    for sample in samples:
        r1_file = f"{base_dir}/{sample}/{sample}_CKDL260002347-1A_23752VLT4_L8_1.fq.gz"
        r2_file = f"{base_dir}/{sample}/{sample}_CKDL260002347-1A_23752VLT4_L8_2.fq.gz"
        
        print(f"Analyzing {sample}...")
        r1_results = check_fastq_lengths(r1_file)
        r2_results = check_fastq_lengths(r2_file)
        
        if r1_results and r2_results:
            all_results[sample] = {
                'R1': r1_results,
                'R2': r2_results
            }
        print("")
    
    # Summary table
    if all_results:
        print("SUMMARY TABLE")
        print("-" * 80)
        print(f"{'Sample':<15} {'R1 Mean':<10} {'R2 Mean':<10} {'R1<=60%':<12} {'R2<=60%':<12}")
        print("-" * 80)
        for sample, results in all_results.items():
            r1_mean = results['R1']['mean']
            r2_mean = results['R2']['mean']
            r1_short = results['R1']['short_pct']
            r2_short = results['R2']['short_pct']
            print(f"{sample:<15} {r1_mean:<10.1f} {r2_mean:<10.1f} {r1_short:<12.1f} {r2_short:<12.1f}")
        
        print("")
        print("INTERPRETATION:")
        print("- If most reads are <= 60 bp, this confirms short fragment library")
        print("- Expected for primer dimers or degraded DNA samples")
        print("- Consistent with high merge rates observed in fastp analysis")


if __name__ == "__main__":
    main()