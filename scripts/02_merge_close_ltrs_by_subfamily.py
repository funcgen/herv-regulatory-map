#!/usr/bin/env python3
"""
Merge nearby LTR BED entries of the same subfamily and strand.
Used to correct for RepeatMasker fragmentation.
"""

import sys
import re

# Maximum gap allowed between LTRs to merge (in bp)
ALLOWED_GAP = 100

def parse_bed_line(line):
    """
    Parse a BED line into a dictionary with:
    - chrom, start, end, name, score, strand
    - subfamily: extracted from the 'name' field as the prefix before '_pos_'
    """
    fields = line.strip().split('\t')
    if len(fields) != 6:
        raise ValueError(f"Invalid BED line (expected 6 fields): {line}")

    chrom, start, end, name, score, strand = fields

    # Use regex to extract subfamily: everything before '_pos_'
    match = re.search(r'^(?P<subfamily>.+?)_pos_', name)
    if match:
        subfamily = match.group('subfamily')
    else:
        subfamily = name  # fallback if no match

    return {
        'chrom': chrom,
        'start': int(start),
        'end': int(end),
        'name': name,
        'score': score,
        'strand': strand,
        'subfamily': subfamily
    }

def format_bed_line(entry):
    """
    Format a merged BED entry back to string, assigning a new standardized name.
    """
    new_name = f"{entry['subfamily']}_merged_pos_{entry['chrom']}_{entry['start']+1}_{entry['end']}_strand_{entry['strand']}"
    return f"{entry['chrom']}\t{entry['start']}\t{entry['end']}\t{new_name}\t{entry['score']}\t{entry['strand']}"

def merge_entries(entries):
    """
    Merge overlapping or nearby entries of the same subfamily on same chrom and strand.
    """
    if not entries:
        return []

    entries.sort(key=lambda x: x['start'])
    merged = [entries[0]]

    for current in entries[1:]:
        last = merged[-1]

        if (current['start'] <= last['end'] + ALLOWED_GAP and
            current['chrom'] == last['chrom'] and
            current['strand'] == last['strand'] and
            current['subfamily'] == last['subfamily']):
            
            last['end'] = max(last['end'], current['end'])
        else:
            merged.append(current)

    return merged

def main(input_file, output_file):
    with open(input_file, 'r') as f:
        entries = [parse_bed_line(line) for line in f if line.strip()]

    # Group by (chrom, strand, subfamily)
    grouped = {}
    for entry in entries:
        key = (entry['chrom'], entry['strand'], entry['subfamily'])
        grouped.setdefault(key, []).append(entry)

    merged_entries = []
    for group in grouped.values():
        merged_entries.extend(merge_entries(group))

    with open(output_file, 'w') as f:
        for entry in merged_entries:
            f.write(format_bed_line(entry) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 02_merge_close_ltrs_by_subfamily.py <input_bed_file> <output_bed_file>")
    else:
        main(sys.argv[1], sys.argv[2])
