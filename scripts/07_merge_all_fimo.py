#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import sys

def merge_fimo_files(input_dir: Path, output_file: Path):
    fimo_files = list(input_dir.glob("*/fimo.tsv"))

    if not fimo_files:
        print(f"❌ No fimo.tsv files found in: {input_dir}")
        sys.exit(1)

    all_dfs = []

    print(f"🔍 Found {len(fimo_files)} FIMO files. Merging with progress bar...")

    for file in tqdm(fimo_files, desc="📦 Processing FIMO files"):
        try:
            df = pd.read_csv(file, sep='\t', comment='#')
            subfamily = file.parent.name
            df.insert(0, "subfamily", subfamily)
            all_dfs.append(df)
        except Exception as e:
            print(f"⚠️ Skipping {file} due to error: {e}")

    merged_df = pd.concat(all_dfs, ignore_index=True)
    merged_df.to_csv(output_file, sep='\t', index=False)

    print(f"✅ Done! Merged {len(fimo_files)} FIMO files into: {output_file}")
    print(f"🧬 Total motif matches: {len(merged_df):,}")


def main():
    parser = argparse.ArgumentParser(
        description="Merge all FIMO result files (fimo.tsv) from subfolders into one file."
    )
    parser.add_argument(
        "--input-dir", "-i", type=Path, required=True,
        help="Path to directory containing FIMO output subfolders"
    )
    parser.add_argument(
        "--output", "-o", type=Path, default=Path("merged_fimo.tsv"),
        help="Output TSV file path (default: merged_fimo.tsv)"
    )
    args = parser.parse_args()

    merge_fimo_files(args.input_dir, args.output)


if __name__ == "__main__":
    main()
