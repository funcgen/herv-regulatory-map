#!/usr/bin/env python3
"""
Parse FIMO TSV output and generate:
1. A TSV file sorted by LTR genomic coordinates.
2. A BED file with the genomic coordinates of each motif hit.

Assumes the sequence name format is:
    <LTR_subfamily>_merged_pos_chr<chromosome>_<start>_<end>_strand_<+|->

Usage:
    python parse_fimo_to_sorted_tsv_and_bed.py \
        --input fimo.tsv \
        --output_tsv fimo_sorted.tsv \
        --output_bed fimo_hits.bed
"""

import pandas as pd
import argparse
import re
from tqdm import tqdm

tqdm.pandas() # enables progress bars for pandas apply/progress_apply


def flip_strand(s):
    """Flip a strand symbol: '+' <-> '-', keep '.' unchanged."""
    if s == '+':
        return '-'
    if s == '-':
        return '+'
    return s

def parse_coordinates_from_name(name):
    """
    Extract chrom, start, end, and strand from sequence names like:
    - LTR2C_merged_pos_chr12_125338183_125338258_strand_-
    - LTR13_pos_GL000224.1_77344_78351
    """
    pattern = r'_(?P<chrom>(chr)?[^_]+)_(?P<start>\d+)_(?P<end>\d+)(?:_strand_(?P<strand>[+-]))?'
    match = re.search(pattern, name)
    if not match:
        raise ValueError(f"Invalid sequence name format: {name}")
    chrom = match.group('chrom')
    start = int(match.group('start'))
    end = int(match.group('end'))
    strand = match.group('strand') if match.group('strand') else '.'
    return chrom, start, end, strand

def main(input_file, output_tsv, output_bed):
    print(f"📥 Reading: {input_file}")
    df = pd.read_csv(input_file, sep='\t', comment='#')

    print(f"🔍 Parsing coordinates from {len(df)} motif hits...")
    coord_data = [parse_coordinates_from_name(name) for name in tqdm(df['sequence_name'], desc="Parsing LTR coordinates")]
    df[['chrom', 'ltr_start', 'ltr_end', 'ltr_strand']] = pd.DataFrame(coord_data, index=df.index)

    print("📊 Sorting TSV by chromosome, LTR start, and motif start...")
    df_sorted = df.sort_values(by=['chrom', 'ltr_start', 'start'])

    print(f"💾 Saving sorted TSV to {output_tsv}")
    df_sorted.to_csv(output_tsv, sep='\t', index=False)

    print("🧬 Generating genomic coordinates for BED...")
    # FIMO start/stop are 1-based inclusive relative to the FASTA sequence.
    # Your FASTA is LTR-oriented (5'->3' for the LTR). For '-' LTRs, the FASTA is reverse-complemented,
    # so we must map using ltr_end (assumed BED end-exclusive).
    plus_mask = (df_sorted['ltr_strand'] == '+')
    minus_mask = (df_sorted['ltr_strand'] == '-')

    # + strand LTR: genomic coordinates increase with FASTA coordinates
    df_sorted.loc[plus_mask, 'genomic_start'] = df_sorted.loc[plus_mask, 'ltr_start'] + (df_sorted.loc[plus_mask, 'start'] - 1)
    df_sorted.loc[plus_mask, 'genomic_end']   = df_sorted.loc[plus_mask, 'ltr_start'] + df_sorted.loc[plus_mask, 'stop']

    # - strand LTR: FASTA is reverse-complemented, so map using ltr_end
    df_sorted.loc[minus_mask, 'genomic_start'] = df_sorted.loc[minus_mask, 'ltr_end'] - df_sorted.loc[minus_mask, 'stop']
    df_sorted.loc[minus_mask, 'genomic_end']   = df_sorted.loc[minus_mask, 'ltr_end'] - (df_sorted.loc[minus_mask, 'start'] - 1)

    print("🧪 Sanity-checking coordinates...")
    bad = df_sorted.progress_apply(lambda r: r["genomic_end"] <= r["genomic_start"], axis=1)
    n_bad = int(bad.sum())
    if n_bad > 0:
        raise ValueError(f"Found {n_bad} hits with genomic_end <= genomic_start. Check coordinate conventions.")

    # BED strand should be genomic: flip FIMO strand when the LTR is on '-'
    df_sorted['genomic_strand'] = df_sorted['strand']
    df_sorted.loc[minus_mask, 'genomic_strand'] = df_sorted.loc[minus_mask, 'strand'].map(flip_strand)


    print("🏷️ Building BED names...")
    df_sorted["bed_name"] = df_sorted.progress_apply(
        lambda r: f"{r['sequence_name']}_{r['motif_id']}_{r['motif_alt_id']}",
        axis=1
    )


    print(f"💾 Saving BED file to {output_bed}")
    bed_df = df_sorted[["chrom", "genomic_start", "genomic_end", "bed_name", "score", "genomic_strand"]]
    bed_df.sort_values(by=["chrom", "genomic_start"]).to_csv(output_bed, sep="\t", header=False, index=False)

    print("✅ Done!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse FIMO TSV and produce sorted TSV and BED.')
    parser.add_argument('--input', required=True, help='Input FIMO TSV file')
    parser.add_argument('--output_tsv', required=True, help='Output sorted TSV file')
    parser.add_argument('--output_bed', required=True, help='Output BED file')
    args = parser.parse_args()
    main(args.input, args.output_tsv, args.output_bed)
