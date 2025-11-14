#!/usr/bin/env python3
"""
Run FIMO to scan LTR sequences for TF binding motifs.
Requires: FIMO from the MEME Suite.

By default, uses raw p-value threshold (not q-values).
"""

import argparse
import subprocess
import os
import sys

def run_fimo(meme_file, fasta_file, output_dir, threshold=1e-4,
             max_scores=1_000_000, use_qval=False, no_pgc=False,
             bgfile=None, fimo_exec="fimo"):

    if not os.path.isfile(meme_file):
        print(f"❌ Motif file not found: {meme_file}")
        sys.exit(1)

    if not os.path.isfile(fasta_file):
        print(f"❌ FASTA file not found: {fasta_file}")
        sys.exit(1)

    if bgfile and not os.path.isfile(bgfile):
        print(f"❌ Background file not found: {bgfile}")
        sys.exit(1)

    if os.path.exists(output_dir):
        print(f"⚠️ Output directory already exists: {output_dir}")
    else:
        os.makedirs(output_dir, exist_ok=True)

    cmd = [
        fimo_exec,
        "--thresh", str(threshold),
        "--max-stored-scores", str(max_scores),
    ]

    if use_qval:
        cmd.append("--qv-thresh")
        print(f"🔬 Using q-value threshold: {threshold}")
    else:
        print(f"🔬 Using raw p-value threshold: {threshold}")

    if no_pgc:
        cmd.append("--no-pgc")

    if bgfile:
        cmd += ["--bgfile", bgfile]
        print(f"🔬 Using background file: {bgfile}")

    cmd += [
        "--oc", output_dir,
        meme_file,
        fasta_file
    ]

    print(f"🚀 Running FIMO on {fasta_file}")
    print(f"🔧 Command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        print(f"✅ FIMO completed. Results in: {output_dir}/fimo.tsv")
    except subprocess.CalledProcessError as e:
        print("❌ FIMO run failed.")
        print(e)
        sys.exit(2)

def main():
    parser = argparse.ArgumentParser(description="Run FIMO on a single FASTA with MEME motif file.")
    parser.add_argument("--motifs", "-m", required=True, help="MEME motif file (e.g. from JASPAR)")
    parser.add_argument("--fasta", "-f", required=True, help="FASTA file with LTR sequences")
    parser.add_argument("--output", "-o", required=True, help="Output directory for FIMO results")
    parser.add_argument("--threshold", "-t", type=float, default=1e-4, help="Threshold (p-value or q-value)")
    parser.add_argument("--qv-thresh", action="store_true", help="Use q-value thresholding instead of p-value")
    parser.add_argument("--no-pgc", action="store_true", help="Suppress parsing coordinates from FASTA headers")
    parser.add_argument("--max-scores", type=int, default=1_000_000, help="Max motif matches to store")
    parser.add_argument("--bgfile", help="Background file generated with fasta-get-markov")
    parser.add_argument("--fimo-path", default="fimo", help="Path to FIMO executable")

    args = parser.parse_args()

    run_fimo(
        meme_file=args.motifs,
        fasta_file=args.fasta,
        output_dir=args.output,
        threshold=args.threshold,
        max_scores=args.max_scores,
        use_qval=args.qv_thresh,
        no_pgc=args.no_pgc,
        bgfile=args.bgfile,
        fimo_exec=args.fimo_path
    )

if __name__ == "__main__":
    main()
