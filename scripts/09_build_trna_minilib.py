#!/usr/bin/env python3
"""
Build a PBS mini-library from GtRNAdb mature tRNA FASTA, with deduplicated outputs.

Inputs
------
- hg38-mature-tRNAs.fa  (from GtRNAdb; mature, intronless tRNAs with CCA)

Outputs
-------
- tRNA_3prime_RNA.fa
- tRNA_3prime_DNA.fa
- tRNA_3prime_DNA_revcomp.fa
- tRNA_3prime_RNA_dedup.fa
- tRNA_3prime_DNA_dedup.fa
- tRNA_3prime_DNA_revcomp_dedup.fa
- manifest.tsv
- dedup_manifest.tsv
- md5sums.txt

Usage
-----
python build_trna_minilib.py \
  --input /path/to/hg38-mature-tRNAs.fa \
  --outdir ./motifs \
  --tail 20
"""
import argparse
import hashlib
from pathlib import Path
import re
import sys
from datetime import datetime

def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def read_fasta(p: Path):
    header = None
    seq_chunks = []
    with p.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if header is not None:
        yield header, "".join(seq_chunks)

def to_dna(rna: str) -> str:
    return rna.upper().replace("U", "T")

def rc_dna(dna: str) -> str:
    tbl = str.maketrans("ACGTN", "TGCAN")
    return dna.translate(tbl)[::-1]

HDR_RE = re.compile(
    r"""
    ^>.*?tRNA-               # ... tRNA-
    (?P<isotype>[A-Za-z]+)   # isotype (e.g., Ala, Lys, Pro)
    -(?P<anticodon>[A-Z]{3}) # anticodon (e.g., AGC)
    """,
    re.X,
)

def parse_hdr(header: str):
    """
    Extract isotype and anticodon if present; fall back to raw header.
    """
    m = HDR_RE.search(header)
    if m:
        return m.group("isotype"), m.group("anticodon")
    return None, None

def main():
    ap = argparse.ArgumentParser(description="Build tRNA PBS mini-library from mature tRNAs.")
    ap.add_argument("--input", required=True, type=Path, help="Path to hg38-mature-tRNAs.fa (GtRNAdb).")
    ap.add_argument("--outdir", required=True, type=Path, help="Output directory (will be created).")
    ap.add_argument("--tail", type=int, default=20, help="Length of 3' tail to extract (default: 20).")
    args = ap.parse_args()

    inp = args.input
    outdir = args.outdir
    tail = args.tail

    if not inp.exists():
        sys.exit(f"[ERROR] Input FASTA not found: {inp}")

    outdir.mkdir(parents=True, exist_ok=True)

    rna_out = outdir / "tRNA_3prime_RNA.fa"
    dna_out = outdir / "tRNA_3prime_DNA.fa"
    rc_out  = outdir / "tRNA_3prime_DNA_revcomp.fa"
    rna_dedup_out = outdir / "tRNA_3prime_RNA_dedup.fa"
    dna_dedup_out = outdir / "tRNA_3prime_DNA_dedup.fa"
    rc_dedup_out  = outdir / "tRNA_3prime_DNA_revcomp_dedup.fa"
    manifest = outdir / "manifest.tsv"
    dedup_manifest = outdir / "dedup_manifest.tsv"
    md5s = outdir / "md5sums.txt"

    n_total = 0
    n_short = 0

    # For deduplication by content
    # Map seq -> info
    rna_uni = {}  # {rna_tail: {"count":int, "ex_header":str, "isotypes":set, "anticodons":set, "full_lens":[], "tail_lens":set}}
    dna_uni = {}
    rc_uni  = {}

    with rna_out.open("w") as frna, dna_out.open("w") as fdna, rc_out.open("w") as frc, manifest.open("w") as fman:
        fman.write("header\tisotype\tanticodon\tfull_len\ttail_len\toutfile_id\n")

        for hdr, seq in read_fasta(inp):
            n_total += 1
            seq = seq.strip().upper()
            if len(seq) < tail:
                n_short += 1
            tail_rna = seq[-tail:] if len(seq) >= tail else seq  # include short ones
            tail_dna = to_dna(tail_rna)
            tail_rc  = rc_dna(tail_dna)

            isotype, anticodon = parse_hdr(hdr)
            out_id = hdr[1:].split()[0]  # safe identifier from header

            # Write raw FASTAs (non-dedup)
            frna.write(f"{hdr}\n{tail_rna}\n")
            fdna.write(f"{hdr}\n{tail_dna}\n")
            frc.write(f"{hdr}_RC\n{tail_rc}\n")

            # Manifest
            fman.write(
                f"{hdr}\t{isotype or ''}\t{anticodon or ''}\t{len(seq)}\t{len(tail_rna)}\t{out_id}\n"
            )

            # Update dedup maps
            def update_uni(d, s):
                if s not in d:
                    d[s] = {
                        "count": 0,
                        "ex_header": hdr,
                        "isotypes": set(),
                        "anticodons": set(),
                        "full_min": len(seq),
                        "full_max": len(seq),
                        "tail_lens": set([len(tail_rna)]),
                    }
                info = d[s]
                info["count"] += 1
                if isotype: info["isotypes"].add(isotype)
                if anticodon: info["anticodons"].add(anticodon)
                info["full_min"] = min(info["full_min"], len(seq))
                info["full_max"] = max(info["full_max"], len(seq))
                info["tail_lens"].add(len(tail_rna))

            update_uni(rna_uni, tail_rna)
            update_uni(dna_uni, tail_dna)
            update_uni(rc_uni,  tail_rc)

    # Write deduplicated FASTAs and their manifest
    with rna_dedup_out.open("w") as frna_d, dna_dedup_out.open("w") as fdna_d, rc_dedup_out.open("w") as frc_d, dedup_manifest.open("w") as fdm:
        fdm.write("space\tsequence\tcount\texample_header\tisotypes\tanticodons\tfull_len_min\tfull_len_max\ttail_lens\n")

        def write_dedup(space_name, uni_map, fh):
            for s in sorted(uni_map.keys()):
                info = uni_map[s]
                # Compact, informative header
                # e.g. >SEQ|count=12|example=tRNA-Lys-CTT|iso={Lys}|anti={CTT}|len={74-76}|tail={19,20}
                iso = ",".join(sorted(info["isotypes"])) if info["isotypes"] else "-"
                anti = ",".join(sorted(info["anticodons"])) if info["anticodons"] else "-"
                lens = ",".join(map(str, sorted(info["tail_lens"])))
                hdr = (
                    f">SEQ|count={info['count']}|example={info['ex_header'][1:].split()[0]}"
                    f"|iso={{{{iso}}}}|anti={{{{anti}}}}|len={info['full_min']}-{info['full_max']}|tail={{{{lens}}}}"
                )
                # fill braces
                hdr = hdr.replace("{{iso}}", iso).replace("{{anti}}", anti).replace("{{lens}}", lens)
                fh.write(f"{hdr}\n{s}\n")

        write_dedup("RNA", rna_uni, frna_d)
        write_dedup("DNA", dna_uni, fdna_d)
        write_dedup("DNA_RC", rc_uni, frc_d)

        # Also dump a table form
        def dump_manifest(space_name, uni_map):
            for s in sorted(uni_map.keys()):
                i = uni_map[s]
                iso = ",".join(sorted(i["isotypes"])) if i["isotypes"] else "-"
                anti = ",".join(sorted(i["anticodons"])) if i["anticodons"] else "-"
                lens = ",".join(map(str, sorted(i["tail_lens"])))
                fdm.write(
                    f"{space_name}\t{s}\t{i['count']}\t{i['ex_header']}\t{iso}\t{anti}\t{i['full_min']}\t{i['full_max']}\t{lens}\n"
                )

        dump_manifest("RNA", rna_uni)
        dump_manifest("DNA", dna_uni)
        dump_manifest("DNA_RC", rc_uni)

    # Checksums (input + outputs) for reproducibility
    with md5s.open("w") as f:
        f.write(f"# Generated: {datetime.utcnow().isoformat()}Z\n")
        for p in [inp, rna_out, dna_out, rc_out, rna_dedup_out, dna_dedup_out, rc_dedup_out, manifest, dedup_manifest]:
            f.write(f"{md5_file(p)}  {p.resolve()}\n")

    print(f"[OK] Processed {n_total} tRNAs (short sequences: {n_short}).")
    print(f"[OK] Wrote:")
    print(f"  - {rna_out}")
    print(f"  - {dna_out}")
    print(f"  - {rc_out}")
    print(f"  - {rna_dedup_out}")
    print(f"  - {dna_dedup_out}")
    print(f"  - {rc_dedup_out}")
    print(f"  - {manifest}")
    print(f"  - {dedup_manifest}")
    print(f"  - {md5s}")

if __name__ == "__main__":
    main()
