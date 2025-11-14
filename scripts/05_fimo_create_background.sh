#!/usr/bin/env bash
# make_ltr_bg.sh — build 0th-order background for FIMO
# Usage:
#   ./make_ltr_bg.sh results/ERV_ltr_v2.fasta bg_LTR_all.txt

set -euo pipefail

IN=${1:-results/ERV_ltr_v2.fasta}
OUT=${2:-bg_LTR_all.txt}

command -v fasta-get-markov >/dev/null 2>&1 || { echo "ERROR: fasta-get-markov not in PATH"; exit 1; }
[ -f "$IN" ] || { echo "ERROR: input FASTA not found: $IN"; exit 1; }

echo "[INFO] Generating 0th-order background from $IN -> $OUT"
fasta-get-markov -m 0 "$IN" "$OUT"

echo "[INFO] Done. Frequencies:"
grep -E '^(A|C|G|T)\s' "$OUT" || true

