#!/usr/bin/env python3
"""
Split an LTR FASTA file into separate FASTAs per subfamily.
Subfamily is defined as the text before the first '_pos_' in the sequence name.

Examples:
  >MLT1K_merged_pos_chr1_21949_22344_strand_+   --> MLT1K
  >LTR12C_pos_chr1_123456_123789_strand_-       --> LTR12C

Output:
  - One FASTA per subfamily
  - Summary TSV with counts

Usage:
  python split_fasta_by_subfamily.py \
    --input input.fasta \
    --output_dir split_by_subfamily
"""

import os
import re
import argparse
from collections import defaultdict
from tqdm import tqdm

def extract_subfamily(header):
    """
    Extracts the subfamily from the FASTA header.
    Everything before '_pos_'. Removes '_merged' suffix if present.
    """
    match = re.search(r'^>(.+?)_pos_', header)
    if match:
        subfam = match.group(1)
        return subfam.replace("_merged", "")
    return header[1:].split()[0]

def parse_fasta_by_subfamily(fasta_path):
    """
    Reads a FASTA file and groups entries by subfamily.
    Returns: dict[subfamily] = list of (header, sequence)
    """
    buckets = defaultdict(list)
    with open(fasta_path) as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and seq_lines:
                    subfam = extract_subfamily(header)
                    buckets[subfam].append((header, ''.join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header and seq_lines:
            subfam = extract_subfamily(header)
            buckets[subfam].append((header, ''.join(seq_lines)))
    return buckets

def write_fastas(buckets, output_dir):
    """
    Writes one FASTA file per subfamily into output_dir.
    Returns a list of (subfamily, count, filename)
    """
    os.makedirs(output_dir, exist_ok=True)
    summary = []
    for subfam, entries in tqdm(buckets.items(), desc="Writing FASTAs"):
        filename = f"{subfam}.fa"
        filepath = os.path.join(output_dir, filename)
        with open(filepath, "w") as out:
            for header, seq in entries:
                out.write(f"{header}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")
        summary.append((subfam, len(entries), filename))
    return summary

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    args = parser.parse_args()

    print(f"📂 Parsing input FASTA: {args.input}")
    buckets = parse_fasta_by_subfamily(args.input)

    print(f"📝 Writing {len(buckets)} subfamily FASTAs to: {args.output_dir}")
    summary = write_fastas(buckets, args.output_dir)

    summary_path = os.path.join(args.output_dir, "subfamily_counts.tsv")
    with open(summary_path, "w") as f:
        f.write("subfamily\tsequence_count\tfile\n")
        for subfam, count, filename in sorted(summary):
            f.write(f"{subfam}\t{count}\t{filename}\n")

    print(f"✅ Done. Summary written to: {summary_path}")

if __name__ == "__main__":
    main()
