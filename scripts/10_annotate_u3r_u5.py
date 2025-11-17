#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
U3/R/U5 annotator (v0.3 - fast)
- Uses pysam.FastaFile (persistent) instead of shelling to samtools per LTR
- Multiprocessing over LTRs
- Precompiled regexes, deduped tRNA tails
- Optional flags to skip PBS/PPT for quick runs
"""

import argparse, sys, re, csv, unicodedata, math, os
from pathlib import Path
from datetime import datetime
from collections import Counter, defaultdict
from bisect import bisect_left, bisect_right
from multiprocessing import Pool, cpu_count
from functools import partial
from itertools import groupby
from tqdm import tqdm

# ---------- Globals for worker processes ----------
FASTA_PATH = None
PYSAM_FASTA = None  # per-process handle
PROMOTER_REGEX = {}
PAS_SET = set()
TRNA_RCS = []
LINK_DIST = 200
U5_MIN_HARD = 40      # reject cleavage if U5 < 40 bp
U5_MAX_SOFT = 320     # ignore very large U5 on sensitive pass
SKIP_PBS = False
SKIP_PPT = False
NAME2CLS = {}
# --- Diagnostics (off by default) ---
DIAG_DPE = False  # emit SIG:DPE_DIAG if any 5-mer DPE-like in +20..+40

# ---- Low-risk promoter tweaks ----
YR_INR_REGEX = r"(?P<YR>[CT][AG])"          # permissive mammalian Inr core at -1/+1 (YR)
YR_INR_WEIGHT = 0.4                         # lighter than your current YYANWYY Inr
GC_BONUS_WINDOW = (-150, +50)               # relative to +1 when scoring a TSS candidate
GC_BONUS_THRESHOLD = 0.60                   # >60% GC counts as "high-GC"
GC_BONUS_NO_TATA = True                     # only apply when no TATA detected near that candidate
SECONDARY_TSS_RADIUS = 25                   # bp around best TSS (for 'broad' mode)
SECONDARY_TSS_FRAC = 0.80                   # accept if score >= 0.8 * best_score



def init_worker(fasta_path, promoter_defs, pas_hex, trna_rc, ctx, link_dist,
                skip_pbs, skip_ppt, name2cls, diag_dpe=False, args_dict=None):
    """Initializer: one FASTA handle per process + compile regexes + motifs."""
    global FASTA_PATH, PYSAM_FASTA, PROMOTER_REGEX, PAS_SET, TRNA_RCS
    global CTX, LINK_DIST, SKIP_PBS, SKIP_PPT, NAME2CLS, DIAG_DPE, args

    # reconstruct lightweight args namespace for workers
    class Args: pass
    args = Args()
    if args_dict:
        for k, v in args_dict.items():
            setattr(args, k, v)
    FASTA_PATH = fasta_path
    try:
        import pysam
    except ImportError as e:
        print("[ERROR] pysam not installed. `conda install -c bioconda pysam`", file=sys.stderr)
        raise
    PYSAM_FASTA = pysam.FastaFile(str(fasta_path))
    PROMOTER_REGEX = {m["name"]: re.compile(m["regex"]) for m in promoter_defs}
    PAS_SET = set(pas_hex)
    # dedupe & keep only >=18 nt
    TRNA_RCS = sorted({s for s in trna_rc if len(s) >= 15}, key=len, reverse=True)
    CTX = ctx
    LINK_DIST = link_dist
    SKIP_PBS = skip_pbs
    SKIP_PPT = skip_ppt
    NAME2CLS = name2cls
    DIAG_DPE = bool(diag_dpe)

# ---------- Basic utils ----------

def rc(seq: str) -> str:
    tbl = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(tbl)[::-1]

def hamming(a, b):
    if len(a) != len(b):
        return math.inf
    return sum(x != y for x, y in zip(a, b))

def read_fasta(path: Path):
    hdr = None; seq = []
    with path.open() as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln: continue
            if ln.startswith(">"):
                if hdr is not None: yield hdr, "".join(seq).upper()
                hdr = ln[1:].split()[0]; seq = []
            else:
                seq.append(ln)
    if hdr is not None:
        yield hdr, "".join(seq).upper()

def load_bed6(path: Path):
    rows=[]
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            chrom, start, end, name, score, strand = parts[:6]
            rows.append({"chrom":chrom,"start":int(start),"end":int(end),
                         "name":name.split()[0],"score":score,"strand":strand})
    return rows

def read_pas_list(path: Path):
    pats=[]
    with path.open() as fh:
        for line in fh:
            s=line.strip().upper()
            if not s or s.startswith("#"): continue
            if not re.fullmatch(r"[ACGT]{6}", s):
                raise ValueError(f"Bad PAS 6-mer: {s}")
            pats.append(s)
    if not pats: raise ValueError("Empty PAS list")
    return pats

def read_promoter_tsv(path: Path):
    motifs = []
    with path.open() as fh:
        for ln_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"): 
                continue
            parts = re.split(r"\t+", line.strip())
            if len(parts) < 4:
                raise ValueError(f"[promoters:{ln_no}] Expected 4 tab-separated fields; got: {line.rstrip()}")
            name, mtype, cons, regex = [x.strip() for x in parts[:4]]

            regex = re.sub(r"\s+", "", regex)
            if not regex.startswith("(?i)"):
                regex = "(?i)" + regex

            if name == "Inr":
                # Add named capture to first literal 'A'
                regex = regex.replace("A", "(?P<A>A)", 1)

            # Validate early
            try:
                re.compile(regex)
            except re.error as e:
                raise ValueError(f"[promoters:{ln_no}] Bad regex for {name}: {regex} :: {e}")

            motifs.append({"name": name, "type": mtype, "consensus": cons, "regex": regex})
    if not motifs:
        raise ValueError("Promoter motif list is empty after parsing.")
    return motifs



def load_trna_tails_rc(path: Path):
    return [seq for _,seq in read_fasta(path)]

# Generic priors (can specialize by ClassFamily if you like)
R_LEN_IDEAL = (60, 350)   # typical R span
U3_MIN_SOFT = 70          # discourage U3 << 70 unless evidence is very strong
R_BELL_WIDTH = 1.15       # widen tolerance by ~15%
PROM_STRONG = 1.4         # threshold for "strong promoter" (Inr ± TATA)

def pick_best_tss(u3_cands, pas_call, L_total, cls_hint=None):
    """
    Select the best U3 end (i.e., TSS-1) by coupling promoter candidates to PAS/cleavage.
    Gentle tweaks:
      - R-length prior slightly widened via R_BELL_WIDTH (lower over-penalization).
      - Ultra-short U3 penalty reduced when promoter evidence is strong (PROM_STRONG).
    """
    if not u3_cands:
        return None, "no_promoter_signal", 0
    if not pas_call or pas_call.get("cleave_idx") is None:
        best = max(u3_cands, key=lambda d: d["score"])
        return best["tss"] - 1, best["evid"], best["score"]

    cleave = pas_call["cleave_idx"]
    best_tuple = None
    mid_R = (R_LEN_IDEAL[0] + R_LEN_IDEAL[1]) / 2.0

    for d in u3_cands:
        tss = d["tss"]
        u3_end_rel = tss - 1
        r_len = cleave - tss

        # R-length prior
        if R_LEN_IDEAL[0] <= r_len <= R_LEN_IDEAL[1]:
            r_ok = 1.0
        else:
            r_ok = max(0.0, 1.0 - (abs(r_len - mid_R) / (400.0 * R_BELL_WIDTH)))

        # U5 gentle bell (~100 nt)
        u5_len = L_total - cleave
        u5_pref = max(0.0, 1.0 - (abs(u5_len - 100) / 200.0))

        # Penalty for ultra-short U3 (reduced if promoter is strong)
        promoter_score = d.get("score", 0.0)
        if u3_end_rel >= U3_MIN_SOFT:
            u3_pen = 0.0
        else:
            u3_pen = 0.25 if promoter_score >= PROM_STRONG else 0.40

        joint = promoter_score + r_ok + 0.5 * u5_pref - u3_pen
        tup = (joint, u3_end_rel, f"{d['evid']};Rlen={r_len};U5len={u5_len}")
        if (best_tuple is None) or (tup[0] > best_tuple[0]):
            best_tuple = tup

    if best_tuple and best_tuple[0] > 0:
        return best_tuple[1], best_tuple[2], best_tuple[0]

    # Fallback: best standalone promoter if coupling didn’t help
    d = max(u3_cands, key=lambda x: x["score"])
    return d["tss"] - 1, d["evid"] + ";no_PAS_coupling", d["score"]


def build_ltr_internal_mapping(ltrs, internals, link_dist=200, strategy="nearest"):
    """
    R-equivalent linker:
      - Require same chrom, strand, and ClassFamily
      - 5' LTR is upstream in transcript space; 3' is downstream
      - Pick nearest (or furthest) within link_dist
    Returns: map_rows (list of dicts), plus helper sets of LTR names seen as 5' and 3'
    """
    # Index internals by (chrom, strand, class)
    by_key = defaultdict(list)
    for iv in internals:
        cls = class_for_name(iv["name"]) or "NA"
        by_key[(iv["chrom"], iv["strand"], cls)].append(iv)
    for key in by_key:
        by_key[key].sort(key=lambda r: (r["start"], r["end"]))

    # group LTRs by (chrom, strand, class) too
    ltrs_by_key = defaultdict(list)
    for lt in ltrs:
        cls = class_for_name(lt["name"]) or "NA"
        ltrs_by_key[(lt["chrom"], lt["strand"], cls)].append(lt)
    for key in ltrs_by_key:
        ltrs_by_key[key].sort(key=lambda r: (r["start"], r["end"]))

    map_rows = []
    names5, names3 = set(), set()

    def pick(cands, which):
        if not cands:
            return None
        if strategy == "furthest":
            cands.sort(key=lambda x: (-x[0], x[1]["name"]))  # largest distance
        else:
            cands.sort(key=lambda x: (x[0], x[1]["name"]))   # NEAREST (R default)
        return cands[0]

    for key, ivs in by_key.items():
        chrom, strand, cls = key
        ltrs_here = ltrs_by_key.get(key, [])
        if not ltrs_here:
            continue

        # Precompute arrays for quick neighbor checks
        l_starts = [r["start"] for r in ltrs_here]
        l_ends   = [r["end"]   for r in ltrs_here]

        for iv in ivs:
            # distances in GENOMIC space, but interpreted in transcript orientation
            # '+' strand: 5' LTR is left (iv.start - ltr.end); 3' LTR is right (ltr.start - iv.end)
            # '-' strand: directions swap
            cand5, cand3 = [], []

            if strand == "+":
                # 5' candidates: LTRs ending before iv.start
                # distance d5 = iv.start - ltr.end
                idx = bisect_left(l_ends, iv["start"])
                lo = max(0, idx-50); hi = min(len(ltrs_here), idx+1)
                for j in range(lo, hi):
                    lt = ltrs_here[j]
                    d5 = iv["start"] - lt["end"]
                    if 0 <= d5 <= link_dist:
                        cand5.append((d5, lt))
                # 3' candidates: LTRs starting after iv.end
                idx = bisect_right(l_starts, iv["end"])
                lo = max(0, idx-1); hi = min(len(ltrs_here), idx+50)
                for j in range(lo, hi):
                    lt = ltrs_here[j]
                    d3 = lt["start"] - iv["end"]
                    if 0 <= d3 <= link_dist:
                        cand3.append((d3, lt))
            else:
                # '-' strand: swap sides in transcript orientation
                # 5' candidates are to the right: d5 = ltr.start - iv.end
                idx = bisect_right(l_starts, iv["end"])
                lo = max(0, idx-1); hi = min(len(ltrs_here), idx+50)
                for j in range(lo, hi):
                    lt = ltrs_here[j]
                    d5 = lt["start"] - iv["end"]
                    if 0 <= d5 <= link_dist:
                        cand5.append((d5, lt))
                # 3' candidates are to the left: d3 = iv.start - ltr.end
                idx = bisect_left(l_ends, iv["start"])
                lo = max(0, idx-50); hi = min(len(ltrs_here), idx+1)
                for j in range(lo, hi):
                    lt = ltrs_here[j]
                    d3 = iv["start"] - lt["end"]
                    if 0 <= d3 <= link_dist:
                        cand3.append((d3, lt))

            pick5 = pick(cand5, "5")
            pick3 = pick(cand3, "3")

            row = {
                "chrom": iv["chrom"],
                "original_start": iv["start"],
                "original_end": iv["end"],
                "strand": iv["strand"],
                "locid": iv["name"],               # internal identifier (R: locid/internal_name)
                "ltr5_name": pick5[1]["name"] if pick5 else ".",
                "dist5":     pick5[0]              if pick5 else None,
                "ltr3_name": pick3[1]["name"] if pick3 else ".",
                "dist3":     pick3[0]              if pick3 else None,
            }
            map_rows.append(row)

            if pick5: names5.add(pick5[1]["name"])
            if pick3: names3.add(pick3[1]["name"])

    return map_rows, names5, names3





# ---------- ClassFamily map ----------

def normalize_repeat_name(s: str) -> str:
    s = s.strip()
    s = unicodedata.normalize("NFKC", s)
    s = s.split()[0]
    if "_pos_" in s: s = s.split("_pos_")[0]
    s = re.sub(r"(_merged)$", "", s, flags=re.IGNORECASE)
    return s

def load_subfamily_map(tsv_path: Path):
    name2cls={}
    with tsv_path.open() as fh:
        rdr=csv.DictReader(fh, delimiter="\t")
        assert "RepeatName" in rdr.fieldnames and "ClassFamily" in rdr.fieldnames, "TSV needs RepeatName,ClassFamily"
        for row in rdr:
            rn=normalize_repeat_name(row["RepeatName"])
            top=row["ClassFamily"].split("/")[-1].upper()
            name2cls[rn]=top
    return name2cls

def class_for_name(name: str):
    return NAME2CLS.get(normalize_repeat_name(name), None)

# ---------- FASTA fetch via pysam ----------

def fetch_seq(chrom: str, start0: int, end0: int) -> str:
    if end0 < start0:
        print(f"[WARN] fetch_seq: end<start for {chrom}:{start0}-{end0}", file=sys.stderr)
    start0 = max(0, start0)
    return PYSAM_FASTA.fetch(chrom, start0, max(start0, end0)).upper()


# ---------- Scanners (precompiled regex used from global) ----------

def scan_promoter_all(seq_plus: str):
    """
    Return TSS candidates with spacing-aware scores using Inr (+ optional TATA),
    plus modest boosts from BREu/BREd (if TATA), DPE (strict +28..+32), MTE (+18..+27),
    SP1 (−150..+10), XCPE1 (−8..+2), and tiny DCE nudges (+6..+34).
    Includes permissive YR Inr fallback and weak TATA-only fallback.
    """
    import re
    L = len(seq_plus)

    # fetch precompiled regexes
    inr_pat   = PROMOTER_REGEX.get("Inr")
    tata_pat  = PROMOTER_REGEX.get("TATA")
    breu_pat  = PROMOTER_REGEX.get("BREu")
    bred_pat  = PROMOTER_REGEX.get("BREd")
    dpe_pat   = PROMOTER_REGEX.get("DPE")
    mte_pat   = PROMOTER_REGEX.get("MTE")
    sp1_pat   = PROMOTER_REGEX.get("SP1")
    xcpe1_pat = PROMOTER_REGEX.get("XCPE1")
    dce1_pat  = PROMOTER_REGEX.get("DCE_SI")
    dce2_pat  = PROMOTER_REGEX.get("DCE_SII")
    dce3_pat  = PROMOTER_REGEX.get("DCE_SIII")

    # weights
    W_TATA_PEAK = 0.80
    W_BREU, W_BRED = 0.15, 0.10
    W_DPE, W_MTE = 0.40, 0.20
    W_SP1, W_XCPE1 = 0.15, 0.10
    W_DCE_STEP = 0.05

    # geometry
    TATA_UP_MIN, TATA_UP_MAX = 10, 60
    TATA_IDEAL, TATA_BW = 30, 20

    def _win(a, b): 
        return max(0, a), min(L, b)

    # helpers for GC bonus
    def _gc_fraction(s):
        if not s:
            return 0.0
        g = s.count("G") + s.count("g")
        c = s.count("C") + s.count("c")
        return (g + c) / float(len(s))

    def _slice_rel(seq, center, rel_lo, rel_hi):
        lo = max(0, center + rel_lo)
        hi = min(len(seq), center + rel_hi)
        return seq[lo:hi]

    # ---------- 1) collect candidate Inr hits ----------
    inr_hits = []
    if inr_pat:
        for m in inr_pat.finditer(seq_plus):
            # Anchor TSS at Inr 'A' if named; else m.start()+2 (YYANWYY)
            try:
                tss = m.start("A")
            except Exception:
                tss = m.start() + 2
            tss = max(1, min(L - 2, tss))
            inr_hits.append({"tss": tss, "label": "Inr_strict", "span": (m.start(), m.end())})

    # ---------- 2) permissive YR Inr fallback ----------
    yr_pat = re.compile(YR_INR_REGEX)
    for m in yr_pat.finditer(seq_plus):
        tss = m.start() + 1  # second base of YR
        inr_hits.append({"tss": tss, "label": "Inr_mam", "span": (m.start(), m.end()),
                         "weight": YR_INR_WEIGHT})

    # deduplicate by TSS position (prefer strict over YR at same site)
    tmp = {}
    for h in inr_hits:
        if h["tss"] not in tmp or h["label"] == "Inr_strict":
            tmp[h["tss"]] = h
    inr_hits = list(tmp.values())

    cands = []

    # ---------- 3) score each Inr-based candidate ----------
    for h in inr_hits:
        tss = h["tss"]
        base = 1.0 * h.get("weight", 1.0)
        
        # Start a fresh evidence list for this TSS candidate
        evid = [h.get("label", "Inr?")]
        
        # Mild penalty for very short U3
        if tss < 60:
            base -= 0.5

        # TATA (position-specific bonus; record span for BRE checks)
        tata_bonus, tata_span = 0.0, None
        if tata_pat:
            up_a, up_b = _win(tss - TATA_UP_MAX, tss - TATA_UP_MIN)
            best = None
            for tm in tata_pat.finditer(seq_plus[up_a:up_b]):
                # distance from TATA end to +1
                dist = tss - (up_a + tm.end())
                sc = max(0.0, 1.0 - abs(dist - TATA_IDEAL) / TATA_BW)
                if best is None or sc > best[0]:
                    best = (sc, dist, (up_a + tm.start(), up_a + tm.end()))
            if best:
                tata_bonus = W_TATA_PEAK * best[0]
                tata_span = best[2]
                evid.append(f"TATAΔ={best[1]}")

        # Downstream elements & BRE
        aux = 0.0

        # DPE: exact 5-nt motif at +28..+32
        if dpe_pat and (tss + 33) <= L:
            if dpe_pat.fullmatch(seq_plus[tss + 28:tss + 33]):
                aux += W_DPE; evid.append("DPE(+28..+32)")

        # MTE: any hit within +18..+27
        if mte_pat and mte_pat.search(seq_plus[tss + 18:tss + 28]):
            aux += W_MTE; evid.append("MTE(+18..+27)")

        # BRE only if TATA is present
        if tata_span:
            ta, tb = tata_span
            if breu_pat and breu_pat.search(seq_plus[max(0, ta - 16):ta]):
                aux += W_BREU; evid.append("BREu")
            if bred_pat and bred_pat.search(seq_plus[tb:min(L, tb + 10)]):
                aux += W_BRED; evid.append("BREd")

        # DCE nudges
        hits = 0
        if dce1_pat and dce1_pat.search(seq_plus[tss + 6:tss + 12]):  hits += 1
        if dce2_pat and dce2_pat.search(seq_plus[tss + 16:tss + 22]): hits += 1
        if dce3_pat and dce3_pat.search(seq_plus[tss + 30:tss + 35]): hits += 1
        if hits:
            aux += W_DCE_STEP * hits; evid.append(f"DCE({hits})")

        # CpG-type helpers
        if sp1_pat and sp1_pat.search(seq_plus[max(0, tss - 150):tss + 10]):
            aux += W_SP1; evid.append("SP1")
        if xcpe1_pat and xcpe1_pat.search(seq_plus[max(0, tss - 8):tss + 2]):
            aux += W_XCPE1; evid.append("XCPE1")

        # Small synergy bonus: SP1 without TATA
        if ("SP1" in evid) and (not tata_span):
            base += 0.05
            evid.append("SP1_noTATA:+0.05")


        # ---- GC% bonus for TATA-less promoters ----
        # expects args.gc_bonus, GC_BONUS_WINDOW, GC_BONUS_THRESHOLD, GC_BONUS_NO_TATA to be defined
        if (globals().get("args") is not None) and getattr(args, "gc_bonus", 0) > 0:
            if (not bool(tata_span)) if GC_BONUS_NO_TATA else True:
                ctx = _slice_rel(seq_plus, tss, GC_BONUS_WINDOW[0], GC_BONUS_WINDOW[1])
                gc = _gc_fraction(ctx)
                if gc >= GC_BONUS_THRESHOLD:
                    base += float(args.gc_bonus)
                    evid.append(f"GCbonus:{gc:.2f}")

        # final score and record
        cands.append({
            "tss": tss,
            "score": base + tata_bonus + aux,
            "has_tata": bool(tata_span),
            "evid": ";".join(evid)
        })

    # ---------- 4) TATA-only fallback ----------
    if not cands and tata_pat:
        for tm in tata_pat.finditer(seq_plus):
            tss = min(L - 2, tm.end() + 30)
            cands.append({
                "tss": tss,
                "score": 0.6,
                "has_tata": True,
                "evid": f"TATA_only@{tm.start()}-{tm.end()}~{tss}"
            })

    return cands






def scan_polyA(seq_plus: str):
    """
    Find a 3′-proximal PAS → cleavage site with optional downstream GU-rich region.
    Strategy:
      - Restrict scan to the last ~220 nt of the LTR.
      - Among candidates, prefer the MOST DOWNSTREAM site that has a clear CA 10–30 nt later.
      - Bonus if [GT]{8,} appears within ~30 nt after the cleavage.
    Output:
      {'cleave_idx': int, 'pas_idx': int, 'score': float, 'evid': str} or None
    """
    L = len(seq_plus)
    if L < 20:
        return None

    WIN = 280
    start = max(0, L - WIN)
    best = None

    # Iterate from 3′ to 5′ to naturally prefer the most downstream valid site
    for i in range(L - 6, start - 1, -1):
        if i < 0: break
        hex6 = seq_plus[i:i+6]
        if hex6 not in PAS_SET:
            continue

        # CA 10–30 nt downstream
        w = seq_plus[i+10:i+31]
        ca_pos = w.find("CA")
        if ca_pos == -1:
            cand = {"cleave_idx": None, "pas_idx": i, "score": 0.8,
                    "evid": f"PAS@{i};no_clear_CA"}
        else:
            ca_abs = i + 10 + ca_pos  # index of 'C' in 'CA'
            # GU-rich 1–30 nt after cleavage (after A of CA)
            gu = re.search(r"[GT]{8,}", seq_plus[ca_abs+1: min(L, ca_abs+1+30)])
            score = 2.2 + (1.0 if gu else 0.0)
            evid  = f"PAS@{i};CA@{ca_abs}" + (f";GU@{(ca_abs+1)+gu.start()}-{(ca_abs+1)+gu.end()}" if gu else "")
            cand = {"cleave_idx": ca_abs+1, "pas_idx": i, "score": score, "evid": evid}

        # Because we iterate from 3′ to 5′, take the first with cleavage
        if cand["cleave_idx"] is not None:
            return cand

        # Otherwise keep the best "no_clear_CA" as fallback
        if best is None or cand["score"] > best["score"]:
            best = cand

    return best


def scan_polyA_sensitive(seq_plus: str):
    """
    Sensitive 3' search with guardrails:
      - broader window (320 nt)
      - CA distance 8..35
      - require U5 >= U5_MIN_HARD (and <= U5_MAX_SOFT) for any candidate
      - A-rich allowed ONLY if a DSE ([GT]{6,}) is present
      - choose best-scoring candidate instead of returning the first
      - weak "no_clear_CA" hint can't force 'OK' (no cleavage index)
    """
    L = len(seq_plus)
    if L < 20:
        return None

    WIN = 360
    start = max(0, L - WIN)

    best = None           # best valid candidate with cleavage
    best_noca = None      # best "no_clear_CA" (weak hint, no cleavage)

    def is_A_rich(h6):
        return h6.count("A") >= 4

    # preference term to keep U5 in a reasonable range (~100 ± 100), gentle
    def u5_pref(u5_len):
        return max(0.0, 0.6 * (1.0 - (abs(u5_len - 100) / 100.0)))

    for i in range(L - 6, start - 1, -1):
        h6 = seq_plus[i:i+6]
        has_hex = (h6 in PAS_SET)
        a_rich  = (not has_hex) and is_A_rich(h6)
        if not (has_hex or a_rich):
            continue

        # broader CA window 8..35
        w = seq_plus[i+8:i+36]
        ca_pos = w.find("CA")

        if ca_pos == -1:
            # only keep a VERY weak hint for A-rich without CA
            if a_rich:
                cand = {"cleave_idx": None, "pas_idx": i, "score": 0.1,
                        "evid": f"A-rich@{i};no_clear_CA"}
                if (best_noca is None) or (cand["score"] > best_noca["score"]):
                    best_noca = cand
            continue

        ca_abs = i + 8 + ca_pos            # index of 'C' in 'CA'
        cleave = ca_abs + 1                # cleavage after A
        u5_len = L - cleave                # how much U5 would remain

        # guardrails on U5 size
        if (u5_len < U5_MIN_HARD) or (u5_len > U5_MAX_SOFT):
            continue

        # DSE
        dse = re.search(r"[GT]{6,}", seq_plus[cleave+10: min(L, cleave + 35)])

        # For A-rich, require DSE; otherwise skip
        if a_rich and not dse:
            continue

        base = 2.4 if has_hex else 1.6
        score = base + (0.6 if dse else 0.0) + u5_pref(u5_len)
        evid  = ("PAS@" if has_hex else "A-rich@") + f"{i};CA@{ca_abs}"
        if dse:
            evid += f";GU@{cleave + dse.start()}-{cleave + dse.end()}"

        cand = {"cleave_idx": cleave, "pas_idx": i, "score": score, "evid": evid}

        # pick the highest score; on ties, prefer the more downstream cleavage
        if (best is None) or (cand["score"] > best["score"]) or \
           (cand["score"] == best["score"] and cleave > best["cleave_idx"]):
            best = cand

    # Prefer a valid cleavage candidate; otherwise return the weak hint (if any)
    return best if best is not None else best_noca




def find_pbs(ds_plus: str, max_mm=None, min_k=12, max_k=20, boundary_bias=True):
    """
    Variable-length fuzzy PBS finder:
      - Try suffixes of each tRNA RC motif from max_k down to min_k (e.g., 20..12)
      - Adaptive mismatch cap: ceil(0.15*k) unless max_mm supplied (uses max of both)
      - Score = (k - mm) - 0.01*offset_from_boundary (if boundary_bias)
    Returns best hit: (s, e, seqpbs, mm, evid) or None
    """
    if SKIP_PBS or not ds_plus:
        return None
    best = None
    L = len(ds_plus)
    for motif_full in TRNA_RCS:
        Kmax = min(max_k, len(motif_full))
        for k in range(Kmax, min_k-1, -1):
            motif = motif_full[-k:]  # use 3' suffix
            mm_cap = max(2, math.ceil(0.15 * k))
            if max_mm is not None:
                mm_cap = max(mm_cap, max_mm)  # allow caller to increase tolerance
            if k > L:
                continue
            for i in range(0, L - k + 1):
                sub = ds_plus[i:i+k]
                mm = hamming(sub, motif)
                if mm <= mm_cap:
                    # prefer closer to LTR boundary (i small)
                    penalty = (0.01 * i) if boundary_bias else 0.0
                    score = (k - mm) - penalty
                    cand = (score, i, i+k, k, mm, motif)
                    if (best is None) or (cand > best):
                        best = cand
    if best:
        score, s, e, k, mm, mot = best
        return s, e, ds_plus[s:e], mm, f"PBS:{mot};len={k};mm={mm}"
    return None


def find_ppt(upstream_seq: str):
    if SKIP_PPT or not upstream_seq:
        return None
    hits=[m for m in re.finditer(r"[AG]{10,}", upstream_seq)]
    if not hits: return None
    L=len(upstream_seq)
    hits.sort(key=lambda m: abs(L - m.end()))
    m=hits[0]
    return m.start(), m.end(), upstream_seq[m.start():m.end()], "PPT"

# ---------- Linking with ClassFamily & (distance = LINK_DIST bp) ----------

def classify_ltr_role(ltr, iv_index):
    """
    iv_index: dict[(chrom,strand,class)] -> dict with 'starts','ends','rows' arrays
    Returns (role, link_name_or_None, reason)
    """
    chrom, strand = ltr["chrom"], ltr["strand"]
    ltr_cls = class_for_name(ltr["name"]) or "NA"
    key = (chrom, strand, ltr_cls)
    if key not in iv_index:
        return "solo", None, "no_internal_for_class"

    idx = iv_index[key]
    starts, ends, rows = idx["starts"], idx["ends"], idx["rows"]

    cand = []
    # Using genomic relations (coordinates), then map to role by strand
    # LTR to the LEFT of internal: ltr.end <= iv.start
    i = bisect_left(starts, ltr["end"])
    for j in range(max(0, i-5), min(len(rows), i+5)):
        d = rows[j]["start"] - ltr["end"]
        if 0 <= d <= LINK_DIST:
            role = "5prime" if strand == "+" else "3prime"
            cand.append((d, role, rows[j]))

    # LTR to the RIGHT of internal: ltr.start >= iv.end
    k = bisect_right(ends, ltr["start"])
    for j in range(max(0, k-5), min(len(rows), k+5)):
        d = ltr["start"] - rows[j]["end"]
        if 0 <= d <= LINK_DIST:
            role = "3prime" if strand == "+" else "5prime"
            cand.append((d, role, rows[j]))

    if not cand:
        return "solo", None, f"no_match_within_{LINK_DIST}bp"

    cand.sort(key=lambda x: x[0])
    best_d, role, iv = cand[0]
    if len(cand) > 1 and (cand[1][0] - best_d) <= 50:
        return "ambiguous", f"{iv['name']}|{cand[1][2]['name']}", f"class={ltr_cls};tie@{best_d}"

    return role, iv["name"], f"class={ltr_cls};strand={strand};d={best_d}"

def build_internal_index(internals):
    """
    Build per (chrom,strand,class) arrays for fast nearest queries.
    """
    buckets=defaultdict(list)
    for iv in internals:
        c = class_for_name(iv["name"]) or "NA"
        buckets[(iv["chrom"], iv["strand"], c)].append(iv)
    index={}
    for key, rows in buckets.items():
        rows_sorted = sorted(rows, key=lambda r: (r["start"], r["end"]))
        starts = [r["start"] for r in rows_sorted]
        ends   = [r["end"]   for r in rows_sorted]
        index[key]={"rows":rows_sorted,"starts":starts,"ends":ends}
    return index

# ---------- Worker function ----------

def process_one(r, ltr_seq_dict, genome_fasta_path):
    """Process a single LTR row. Uses global PYSAM_FASTA & regexes set in init()."""
    # --- unpack row ---
    name   = r["name"]
    chrom  = r["chrom"]
    strand = r["strand"]
    start  = r["start"]
    end    = r["end"]

    seq = ltr_seq_dict.get(name, "")
    if not seq:
        return None  # no sequence available; skip

    role        = r.get("role", "solo")
    link_name   = r.get("link_internal", None)
    reason      = r.get("link_reason", "")

    # -------------------------------------------------------------------------
    # LTR body is already strand-corrected in your FASTA generator.
    # IMPORTANT: do NOT RC here again; keep a single transcript (+) orientation.
    # -------------------------------------------------------------------------
    seq_plus = seq
    L = len(seq_plus)

    # --- signals INSIDE the LTR (promoter & initial PAS scan) ---
    u3_cands = scan_promoter_all(seq_plus)

    # First pass: permissive scan over the whole LTR (used to help rank TSS)
    pas_call_full = scan_polyA(seq_plus)
    if not pas_call_full or pas_call_full.get("cleave_idx") is None:
        pas_call_full = scan_polyA_sensitive(seq_plus)

    # Pick best TSS using whole-LTR PAS evidence just for scoring
    u3_end_rel, prom_evid, prom_score = pick_best_tss(
        u3_cands, pas_call_full, L, class_for_name(name)
    )

    # -------------------------------------------------------------------------
    # ENFORCE PAS DOWNSTREAM OF TSS IN TRANSCRIPT SPACE
    # Re-scan ONLY the tail after TSS; if none found, accept full-scan PAS
    # only if it lies downstream. Otherwise, treat PAS as absent.
    # -------------------------------------------------------------------------
    pas_call = None
    if u3_end_rel is not None and u3_end_rel >= 0 and u3_end_rel < L - 1:
        tail_off = u3_end_rel + 1
        tail = seq_plus[tail_off:]
        pas_tail = scan_polyA(tail)
        if not pas_tail or pas_tail.get("cleave_idx") is None:
            pas_tail = scan_polyA_sensitive(tail)

        if pas_tail and pas_tail.get("cleave_idx") is not None:
            # Offset the indices to full LTR coordinates
            pas_call = dict(pas_tail)
            pas_call["cleave_idx"] = tail_off + pas_tail["cleave_idx"]
            if "pas_idx" in pas_tail and pas_tail["pas_idx"] is not None:   # <-- FIX
                pas_call["pas_idx"] = tail_off + pas_tail["pas_idx"]        # <-- FIX
            # keep score/evidence; optionally annotate origin
            pas_call["evid"] = f"{pas_tail.get('evid','PAS')}@tail+{tail_off}"
        else:
            # No PAS in tail: use full-scan PAS only if downstream of TSS
            if pas_call_full and pas_call_full.get("cleave_idx") is not None \
               and pas_call_full["cleave_idx"] > u3_end_rel:
                pas_call = pas_call_full
    else:
        # No valid TSS: fall back to best PAS from full scan (may be None)
        pas_call = pas_call_full if (pas_call_full and pas_call_full.get("cleave_idx") is not None) else None

    # Resolve U5 start & PAS fields
    def _strip_tokens(ev: str, keys=("Rlen","U5len")) -> str:
        """Remove any existing Rlen=/U5len= tokens from an evidence string."""
        if not ev:
            return ev
        pat = re.compile(r"(?:^|;)(?:" + "|".join(keys) + r")=[^;]+")
        ev2 = pat.sub("", ev)
        ev2 = re.sub(r";{2,}", ";", ev2).strip(";")
        return ev2

    if pas_call and pas_call.get("cleave_idx") is not None:
        u5_start_rel = pas_call["cleave_idx"]
        pas_evid     = pas_call.get("evid", "PAS")
        pas_score    = pas_call.get("score", 0)
    else:
        u5_start_rel = None
        pas_evid     = "no_PAS_signal"
        pas_score    = 0

    # --- fetch flanks (genomic, then normalize to transcript '+') ---
    if strand == "+":
        ds_start0, ds_end0 = end, end + CTX             # downstream of LTR end (5' LTR PBS context)
        us_start0, us_end0 = max(0, start - CTX), start # upstream of LTR start (3' LTR PPT context)
    else:
        # For '-' LTRs, transcript-downstream is LEFT in genome, upstream is RIGHT
        ds_start0, ds_end0 = max(0, start - CTX), start
        us_start0, us_end0 = end, end + CTX

    ds_seq = fetch_seq(chrom, ds_start0, ds_end0)
    us_seq = fetch_seq(chrom, us_start0, us_end0)

    # Normalize flanks to transcript (+) orientation
    ds_plus = ds_seq if strand == "+" else rc(ds_seq)
    us_plus = us_seq if strand == "+" else rc(us_seq)

    # --- PBS (immediately downstream of 5' LTR) ---
    pbs_call = None
    if role in ("5prime", "ambiguous",  "both") and not SKIP_PBS and ds_plus:
        PBS_MAXWIN = 250  # scan the first ~250 bp
        ds_scan = ds_plus[:PBS_MAXWIN]
        pbs_call = find_pbs(ds_scan, max_mm=None, min_k=12, max_k=20, boundary_bias=True)
        # -> (s, e, seqpbs, mm, evid_str) or None

    # --- PPT (immediately upstream of 3' LTR) ---
    ppt_call = None
    if role in ("3prime", "ambiguous", "both") and not SKIP_PPT and us_plus:
        ppt_call = find_ppt(us_plus)  # -> (s, e, seqppt, "PPT") or None

    # Prepare outputs
    segs_ok = []       # only OK segments (U3/R/U5)
    segs_all = []      # includes LOW_CONF/partial
    flank_segs = []    # PBS/PPT in flanks.bed
    signal_segs = []   # CLEAN signals actually used (this function now)

    # Helper: map transcriptional '+' window indices back to genomic
    def plus_idx_to_genome(win_start0, win_end0, i, j, strand_is_plus):
        """Convert [i,j) within fetched [win_start0,win_end0) back to genomic [g0,g1)."""
        if strand_is_plus:
            return (win_start0 + i, win_start0 + j)
        else:
            # reverse-complemented window: [i,j) -> [win_end0-j, win_end0-i)
            return (win_end0 - j, win_end0 - i)

    # Map PBS to genomic coords (FLANKS only; SIGNALS will add one clean line later)
    pbs_sig = None
    if pbs_call:
        s, e, seqpbs, mm, evid_pbs = pbs_call  # indices in ds_plus
        g_start, g_end = plus_idx_to_genome(ds_start0, ds_end0, s, e, strand == "+")
        pbs_score = max(0, (e - s) - mm)
        flank_segs.append((chrom, g_start, g_end, f"{name}|PBS", pbs_score, strand))
        pbs_sig = (g_start, g_end, pbs_score)

    # Map PPT to genomic coords (FLANKS only; SIGNALS will add one clean line later)
    ppt_sig = None
    if ppt_call:
        s, e, seqppt, evid_ppt = ppt_call  # indices in us_plus
        g_start, g_end = plus_idx_to_genome(us_start0, us_end0, s, e, strand == "+")
        ppt_score = (e - s)
        flank_segs.append((chrom, g_start, g_end, f"{name}|PPT", ppt_score, strand))
        ppt_sig = (g_start, g_end, ppt_score)

    # Confidence & status
    conf = prom_score + pas_score + (1 if pbs_call else 0) + (1 if ppt_call else 0)
    status = "OK"
    if (u3_end_rel is None) or (u5_start_rel is None) or (u3_end_rel >= u5_start_rel) or (u3_end_rel < 0) or (u5_start_rel > L):
        status = "LOW_CONF"

    # Evidence string
    # Clean promoter evidence from stale geometry tokens
    prom_evid_clean = _strip_tokens(prom_evid)

    # Recompute FINAL geometry tokens from FINAL boundaries
    final_tokens = []
    # u3_end_rel and u5_start_rel must already be resolved at this point
    if (u3_end_rel is not None) and (u5_start_rel is not None) and (u5_start_rel > u3_end_rel):
        # Match your R analysis convention: R_len = u5_start_rel - u3_end_rel
        Rlen_final = int(u5_start_rel - u3_end_rel)
        # U5 extends from u5_start_rel to the end of the LTR
        # 'L' should be the LTR length already defined in your code (len(seq_plus))
        U5len_final = int(L - u5_start_rel)
        final_tokens.append(f"Rlen={Rlen_final}")
        final_tokens.append(f"U5len={U5len_final}")

    # Build the final evidence string in this order:
    #  1) promoter evidence (cleaned)
    #  2) final geometry tokens (Rlen/U5len) if available
    #  3) PAS evidence (from the downstream-of-TSS-checked pass)
    #  4) PBS/PPT extras
    evid_parts = [prom_evid_clean]
    if final_tokens:
        evid_parts.append(";".join(final_tokens))
    evid_parts.append(pas_evid)
    if pbs_call:
        evid_parts.append(f"{pbs_call[4]}@ds[{pbs_call[0]}-{pbs_call[1]}]")
    if ppt_call:
        evid_parts.append(f"{ppt_call[3]}@us[{ppt_call[0]}-{ppt_call[1]}]")

    evidence = ";".join([p for p in evid_parts if p])


    # Map plus-oriented spans back to genome
    def plus_span_to_genome(a, b):
        if strand == "+":
            return (start + a, start + b)
        else:
            # map [a, b) in plus to genomic [end-b, end-a)
            return (end - b, end - a)

    # Build U3/R/U5 segments
    if status == "OK":
        # plus-relative half-open spans
        U3_span = (0, u3_end_rel + 1)
        R_span  = (u3_end_rel + 1, u5_start_rel)
        U5_span = (u5_start_rel, L)

        U3_start0, U3_end0 = plus_span_to_genome(*U3_span)
        R_start0,  R_end0  = plus_span_to_genome(*R_span)
        U5_start0, U5_end0 = plus_span_to_genome(*U5_span)

        segs_ok = [
            (chrom, U3_start0, U3_end0, f"{name}|U3", conf, strand),
            (chrom, R_start0,  R_end0,  f"{name}|R",  conf, strand),
            (chrom, U5_start0, U5_end0, f"{name}|U5", conf, strand),
        ]
        segs_all.extend(segs_ok)
    else:
        # Build partial segments when boundaries are incomplete
        plus_spans = []
        if (u3_end_rel is not None) and (u3_end_rel >= 0):
            plus_spans.append(("U3", 0, u3_end_rel + 1))
        if (u3_end_rel is not None) and (u5_start_rel is not None) and (u5_start_rel > u3_end_rel + 1):
            plus_spans.append(("R", u3_end_rel + 1, u5_start_rel))
        if (u5_start_rel is not None) and (u5_start_rel < L):
            plus_spans.append(("U5", u5_start_rel, L))

        for lab, a, b in plus_spans:
            g0, g1 = plus_span_to_genome(a, b)
            segs_all.append((chrom, g0, g1, f"{name}|{lab}|LOW_CONF", conf, strand))


    # =========================
    # SYMMETRIC "USED-ONLY" SIGNALS
    # =========================
    def emit_symmetric_signals():
        """Emit symmetric, used signals for ALL LTRs (5′, 3′, solo), if present."""
        best_tss = (u3_end_rel + 1) if (u3_end_rel is not None and u3_end_rel >= 0) else None

        def to_genomic(a, b):
            if strand == "+":
                return (start + a, start + b)
            else:
                return (end - b, end - a)

        # ---- Promoter/TSS block (inside LTR)
        if best_tss is not None:
            # Explicit TSS point
            g0, g1 = to_genomic(best_tss, best_tss + 1)
            signal_segs.append((chrom, g0, g1, f"{name}|SIG:TSS", 1, strand))
            # Secondary TSS diagnostics (broad mode only) — does NOT affect U3/R/U5 segmentation
            if getattr(args, "tss_mode", "focused") == "broad" and best_tss is not None and u3_cands:
                sec = []
                for c in u3_cands:
                    # skip the primary
                    if c.get("tss") == best_tss:
                        continue
                    # within radius and strong enough relative to best
                    if abs(c["tss"] - best_tss) <= SECONDARY_TSS_RADIUS and c["score"] >= SECONDARY_TSS_FRAC * float(prom_score):
                        # derive a short label from evidence if needed
                        lab = c.get("label")
                        if not lab:
                            lab = (c.get("evid", "").split(";")[0] or "Inr?")
                        sec.append((c, lab))

                # emit as diagnostic/auxiliary signals
                for c, lab in sec:
                    a = int(c["tss"]); b = a + 1
                    g0, g1 = to_genomic(a, b)
                    signal_segs.append((
                        chrom, g0, g1,
                        f"{name}|SIG:TSS_SECONDARY|{lab}",
                        round(float(c.get("score", 0.0)), 3),
                        strand
                    ))


            inr_pat   = PROMOTER_REGEX.get("Inr")
            tata_pat  = PROMOTER_REGEX.get("TATA")
            breu_pat  = PROMOTER_REGEX.get("BREu")
            bred_pat  = PROMOTER_REGEX.get("BREd")
            dpe_pat   = PROMOTER_REGEX.get("DPE")
            mte_pat   = PROMOTER_REGEX.get("MTE")
            sp1_pat   = PROMOTER_REGEX.get("SP1")
            xcpe1_pat = PROMOTER_REGEX.get("XCPE1")

            # Inr: prefer the match whose A+1 equals best_tss; otherwise nearest A to best_tss
            if inr_pat:
                w0 = max(0, best_tss - 6); w1 = min(L, best_tss + 6)
                best = None
                for mm in inr_pat.finditer(seq_plus[w0:w1]):
                    a = w0 + mm.start()
                    b = w0 + mm.end()
                    # try named A group; else fallback to first 'A' position (+2)
                    if mm.re.groupindex.get('A'):
                        Apos = w0 + mm.start('A')
                    else:
                        Apos = a + 2
                    # exact A alignment wins immediately
                    if Apos == best_tss:
                        best = ((0, 0), a, b); break
                    key = (abs(Apos - best_tss), -(b - a))
                    if (best is None) or (key < best[0]):
                        best = (key, a, b)
                if best:
                    _, a, b = best
                    g0, g1 = to_genomic(a, b)
                    signal_segs.append((chrom, g0, g1, f"{name}|SIG:Inr", 1, strand))

            # TATA: 10–60 upstream (closest to −30)
            tata_span = None
            if tata_pat:
                u0 = max(0, best_tss - 60); u1 = max(0, best_tss - 10)
                best = None
                for m in tata_pat.finditer(seq_plus[u0:u1]):
                    a = u0 + m.start(); b = u0 + m.end()
                    dist = best_tss - b   # distance from TATA end to +1
                    if 10 <= dist <= 60:
                        key = (abs(dist - 30),)
                        if (best is None) or (key < best[0]):
                            best = (key, a, b, dist)
                if best:
                    _, a, b, _ = best
                    tata_span = (a, b)
                    g0, g1 = to_genomic(a, b)
                    signal_segs.append((chrom, g0, g1, f"{name}|SIG:TATA", 1, strand))

            # BREu (tight: [TATA_start−16, TATA_start)) and BREd ([TATA_end, TATA_end+10))
            if tata_span:
                ta, tb = tata_span
                if breu_pat:
                    u0 = max(0, ta - 16); u1 = ta
                    best = None
                    for m in breu_pat.finditer(seq_plus[u0:u1]):
                        a = u0 + m.start(); b = u0 + m.end()
                        dist = ta - a
                        if 0 < dist <= 16:
                            key = (abs(dist - 8),)
                            if (best is None) or (key < best[0]):
                                best = (key, a, b)
                    if best:
                        _, a, b = best
                        g0, g1 = to_genomic(a, b)
                        signal_segs.append((chrom, g0, g1, f"{name}|SIG:BREu", 1, strand))
                if bred_pat:
                    d0 = tb; d1 = min(L, tb + 10)
                    best = None
                    for m in bred_pat.finditer(seq_plus[d0:d1]):
                        a = d0 + m.start(); b = d0 + m.end()
                        dist = a - tb
                        if 1 <= dist <= 7:
                            key = (abs(dist - 4),)
                            if (best is None) or (key < best[0]):
                                best = (key, a, b)
                    if best:
                        _, a, b = best
                        g0, g1 = to_genomic(a, b)
                        signal_segs.append((chrom, g0, g1, f"{name}|SIG:BREd", 1, strand))

            # DPE: strict — fullmatch 5 nt at +28..+32 from A+1
            if dpe_pat and (best_tss + 33) <= L:
                a = best_tss + 28; b = best_tss + 33
                if dpe_pat.fullmatch(seq_plus[a:b]):
                    g0, g1 = to_genomic(a, b)
                    signal_segs.append((chrom, g0, g1, f"{name}|SIG:DPE", 1, strand))
            
            # DPE_DIAG: optional diagnostic — best 5-mer anywhere in +20..+40 (no scoring impact)
            if DIAG_DPE and dpe_pat:
                # scan +20..+40, pick the hit whose center is closest to +30
                w0 = min(L, best_tss + 20); w1 = min(L, best_tss + 40)
                best = None
                for m in dpe_pat.finditer(seq_plus[w0:w1]):
                    a = w0 + m.start(); b = w0 + m.end()
                    center = (a + b) / 2.0
                    key = (abs((best_tss + 30) - center), -(b - a))  # prefer ~+30, then longer (defensive)
                    if (best is None) or (key < best[0]):
                        best = (key, a, b)
                if best:
                    _, a, b = best
                    g0, g1 = to_genomic(a, b)
                    signal_segs.append((chrom, g0, g1, f"{name}|SIG:DPE_DIAG", 1, strand))

            # SP1 / GC-box: one best in [-150..+10] (prefer around −30)
            if sp1_pat:
                u0 = max(0, best_tss - 150); u1 = min(L, best_tss + 10)
                best = None
                for m in sp1_pat.finditer(seq_plus[u0:u1]):
                    a = u0 + m.start(); b = u0 + m.end()
                    center = (a + b) / 2.0
                    key = (abs((best_tss - 30) - center), -(b - a))
                    if (best is None) or (key < best[0]):
                        best = (key, a, b)
                if best:
                    _, a, b = best
                    g0, g1 = to_genomic(a, b)
                    signal_segs.append((chrom, g0, g1, f"{name}|SIG:SP1", 1, strand))

            # MTE: any occurrence in +18..+27; report the closest-to-+22 center
            if mte_pat:
                m0 = min(L, best_tss + 18); m1 = min(L, best_tss + 27)
                best = None
                for m in mte_pat.finditer(seq_plus[m0:m1]):
                    a = m0 + m.start(); b = m0 + m.end()
                    key = (abs((best_tss + 22) - ((a + b) / 2.0)), -(b - a))
                    if (best is None) or (key < best[0]): best = (key, a, b)
                if best:
                    _, a, b = best
                    g0, g1 = to_genomic(a, b)
                    signal_segs.append((chrom, g0, g1, f"{name}|SIG:MTE", 1, strand))

            # DCE subelements: SI (+6..+11), SII (+16..+21), SIII (+28..+34) — one per subelement
            for lab, pat, rel0, rel1 in (
                ("DCE_SI",  PROMOTER_REGEX.get("DCE_SI"),   6, 11),
                ("DCE_SII", PROMOTER_REGEX.get("DCE_SII"), 16, 21),
                ("DCE_SIII",PROMOTER_REGEX.get("DCE_SIII"),28, 34),
            ):
                if pat:
                    w0 = max(0, best_tss + rel0); w1 = min(L, best_tss + rel1 + 1)
                    for m in pat.finditer(seq_plus[w0:w1]):
                        a = w0 + m.start(); b = w0 + m.end()
                        g0, g1 = to_genomic(a, b)
                        signal_segs.append((chrom, g0, g1, f"{name}|SIG:{lab}", 1, strand))
                        break  # one is enough per subelement

            # XCPE1 (−8..+2) — report one if present
            if xcpe1_pat:
                x0 = max(0, best_tss - 8); x1 = min(L, best_tss + 2)
                for m in xcpe1_pat.finditer(seq_plus[x0:x1]):
                    a = x0 + m.start(); b = x0 + m.end()
                    g0, g1 = to_genomic(a, b)
                    signal_segs.append((chrom, g0, g1, f"{name}|SIG:XCPE1", 1, strand))
                    break

        # ---- PAS + CLEAVAGE (inside LTR)
        if pas_call:
            i = pas_call.get("pas_idx")
            if i is not None:
                g0, g1 = to_genomic(i, i + 6)
                signal_segs.append((chrom, g0, g1, f"{name}|SIG:PAS", int(pas_call.get("score", 1)), strand))
            cleave = pas_call.get("cleave_idx")
            if cleave is not None:
                g0, g1 = to_genomic(cleave, cleave + 1)
                signal_segs.append((chrom, g0, g1, f"{name}|SIG:CLEAVAGE", int(pas_call.get("score", 1)), strand))

        # ---- Flanks (PBS/PPT) — emit if found (no role gating)
        if pbs_sig:
            g0, g1, sc = pbs_sig
            signal_segs.append((chrom, g0, g1, f"{name}|SIG:PBS", sc, strand))
        if ppt_sig:
            g0, g1, sc = ppt_sig
            signal_segs.append((chrom, g0, g1, f"{name}|SIG:PPT", sc, strand))



    # Emit canonical signals to all LTRs
    emit_symmetric_signals()

    row = {
        "name": name, "chrom": chrom, "start": start, "end": end, "strand": strand,
        "role": role, "link_internal": (link_name or ""), "link_reason": reason,
        "u3_end_rel": ("" if u3_end_rel is None else u3_end_rel),
        "u5_start_rel": ("" if u5_start_rel is None else u5_start_rel),
        "confidence": conf, "status": status, "evidence": evidence
    }
    return row, segs_ok, segs_all, flank_segs, signal_segs



# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(description="Annotate U3/R/U5 in LTRs (fast, ClassFamily-aware).")
    ap.add_argument("--ltr-bed", required=True, type=Path)
    ap.add_argument("--ltr-fasta", required=True, type=Path)
    ap.add_argument("--genome", required=True, type=Path)
    ap.add_argument("--pas", required=True, type=Path)
    ap.add_argument("--promoters", required=True, type=Path)
    ap.add_argument("--tss-mode",
                choices=("focused", "broad"),
                default="focused",
                help="Focused (single +1) or broad (emit nearby secondary TSS signals). Default: focused.")
    ap.add_argument("--gc-bonus",
                    type=float,
                    default=0.15,
                    help="Small score bonus (0–0.5) for high-GC TATA-less promoters. 0 disables.")
    ap.add_argument("--trna-rc", required=True, type=Path)
    ap.add_argument("--internal-bed", required=True, type=Path)
    ap.add_argument("--subfamily-map", required=True, type=Path)
    ap.add_argument("--link-dist", type=int, default=200) 
    ap.add_argument("--link-strategy", choices=("nearest","furthest"), default="nearest",
                help="When multiple candidates exist within link-dist, pick the nearest (R default) or furthest.")
    ap.add_argument("--ctx", type=int, default=200)
    ap.add_argument("--threads", type=int, default=max(1, cpu_count()//2))
    ap.add_argument("--no-pbs", action="store_true", help="Skip PBS search (faster)")
    ap.add_argument("--no-ppt", action="store_true", help="Skip PPT search (faster)")
    ap.add_argument("--outdir", type=Path, default=Path("outputs"))
    ap.add_argument("--diag-dpe", action="store_true",
                    help="Emit SIG:DPE_DIAG if any DPE-like 5-mer occurs in +20..+40")
    args = ap.parse_args()

    outdir=args.outdir; (outdir/"logs").mkdir(parents=True, exist_ok=True)
    logf = outdir/"logs"/f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    log = open(logf,"w"); print("## U3/R/U5 annotator v0.3", file=log)

    # Load resources (in parent proc)
    pas_hex = read_pas_list(args.pas)
    prom_defs = read_promoter_tsv(args.promoters)
    trna_rc = load_trna_tails_rc(args.trna_rc)
    name2cls = load_subfamily_map(args.subfamily_map)

    print(f"PAS hexamers: {len(pas_hex)}", file=log)
    print(f"Promoter motifs: {len(prom_defs)}", file=log)
    for m in prom_defs:
        print(f"motif {m['name']}: {m['regex']}", file=log)

    print(f"tRNA RC tails (raw): {len(trna_rc)}", file=log)
    print(f"ClassFamily entries: {len(name2cls)}", file=log)

    ltrs = load_bed6(args.ltr_bed)
    internals = load_bed6(args.internal_bed)
    ltr_fa = {h:s for h,s in read_fasta(args.ltr_fasta)}

    # Build internal index
    global NAME2CLS
    NAME2CLS = name2cls
    iv_index = build_internal_index(internals)

    prepped = [r.copy() for r in ltrs]


    # Build mapping in the parent process (R-equivalent)
    map_rows, names5, names3 = build_ltr_internal_mapping(
        ltrs, internals, link_dist=args.link_dist, strategy=args.link_strategy
    )

    # Assign roles from the mapping sets (exact R behavior):
    #  - if an LTR appears as 5' for any internal → mark 5prime
    #  - if it appears as 3' for any internal → mark 3prime
    #  - if both → "both"
    name2role = {}
    for nm in set([r["name"] for r in ltrs]):
        has5 = nm in names5
        has3 = nm in names3
        if has5 and has3:
            name2role[nm] = "both"
        elif has5:
            name2role[nm] = "5prime"
        elif has3:
            name2role[nm] = "3prime"
        else:
            name2role[nm] = "solo"

    # Attach role & link info to rows we pass to workers
    # choose a single 'link_internal' pointer preferentially from 5' then 3', nearest
    # (R picks per-internal; here we pick a stable pointer per LTR for bookkeeping)
    from collections import defaultdict
    ltr2cands = defaultdict(list)
    for row in map_rows:
        if row["ltr5_name"] != ".":
            ltr2cands[row["ltr5_name"]].append(("5", row["locid"], row["dist5"]))
        if row["ltr3_name"] != ".":
            ltr2cands[row["ltr3_name"]].append(("3", row["locid"], row["dist3"]))

    prepped = []
    for r in ltrs:
        nm = r["name"]
        role = name2role.get(nm, "solo")
        link_name = ""
        link_reason = f"role={role};class={class_for_name(nm) or 'NA'}"

        if nm in ltr2cands:
            # prefer 5' over 3'; then nearest distance
            c = sorted(ltr2cands[nm], key=lambda x: (0 if x[0]=="5" else 1, x[2] if x[2] is not None else 10**9))
            link_name = c[0][1]
            link_reason += f";d={c[0][2]}"

        r2 = r.copy()
        r2["role"] = role
        r2["link_internal"] = link_name
        r2["link_reason"] = link_reason
        prepped.append(r2)




    # Init workers
    init_kwargs = dict(
        fasta_path=args.genome,
        promoter_defs=prom_defs,
        pas_hex=pas_hex,
        trna_rc=trna_rc,
        ctx=args.ctx,
        link_dist=args.link_dist,
        skip_pbs=args.no_pbs,
        skip_ppt=args.no_ppt,
        name2cls=name2cls,
        diag_dpe=args.diag_dpe,
        args_dict={
            "tss_mode": args.tss_mode,
            "gc_bonus": args.gc_bonus
        }
    )

    with Pool(processes=args.threads, initializer=init_worker, initargs=tuple(init_kwargs.values())) as pool:
        worker = partial(process_one, ltr_seq_dict=ltr_fa, genome_fasta_path=args.genome)
        results=[]
        for out in tqdm(pool.imap_unordered(worker, prepped, chunksize=200), total=len(prepped), unit="ltr", desc="Annotating LTRs"):
            if out is not None:
                results.append(out)

    # Collect
    rows = []
    segs_ok = []
    segs_all = []
    flank = []
    signals = []

    for out in results:
        row, seg_ok, seg_all, fseg, sigs = out
        rows.append(row)
        segs_ok.extend(seg_ok)
        segs_all.extend(seg_all)
        flank.extend(fseg)
        signals.extend(sigs)

    # Write outputs
    map_path = outdir / "ERV_full_plus_components.map.tsv"
    bed_components_path = outdir / "ERV_full_plus_components.bed"

    # 1) Mapping table (R: map_dt)
    with map_path.open("w") as f:
        f.write("\t".join(["chrom","original_start","original_end","strand","locid",
                        "ltr5_name","dist5","ltr3_name","dist3"]) + "\n")
        for m in map_rows:
            f.write("\t".join(map(str, [
                m["chrom"], m["original_start"], m["original_end"], m["strand"], m["locid"],
                m["ltr5_name"], ("" if m["dist5"] is None else m["dist5"]),
                m["ltr3_name"], ("" if m["dist3"] is None else m["dist3"])
            ])) + "\n")

    # 2) Components BED (optional, mirrors R helper)
    with bed_components_path.open("w") as f:
        # internal
        for iv in internals:
            f.write(f"{iv['chrom']}\t{iv['start']}\t{iv['end']}\t{iv['name']}\t0\t{iv['strand']}\n")
        # any LTR that appears in mapping (either side)
        linked = set([m["ltr5_name"] for m in map_rows if m["ltr5_name"] != "."] +
                    [m["ltr3_name"] for m in map_rows if m["ltr3_name"] != "."])
        for lt in ltrs:
            if lt["name"] in linked:
                f.write(f"{lt['chrom']}\t{lt['start']}\t{lt['end']}\t{lt['name']}\t0\t{lt['strand']}\n")





    tsv_path = outdir / "U3R_U5_catalogue.tsv"
    bed_ok_path = outdir / "U3R_U5_segments.bed"          # OK-only (as before)
    bed_all_path = outdir / "U3R_U5_segments_all.bed"     # includes LOW_CONF & partial
    flank_path = outdir / "U3R_U5_flanks.bed"             # PBS/PPT near LTRs
    signals_path = outdir / "U3R_U5_signals.bed"          # promoter/PAS/cleavage + (PBS/PPT also added)

    with tsv_path.open("w") as tsv:
        tsv.write("\t".join(["name","chrom","start","end","strand",
                             "role","link_internal","link_reason",
                             "u3_end_rel","u5_start_rel","confidence","status","evidence"])+"\n")
        for r in rows:
            tsv.write("\t".join(map(str, [
                r["name"], r["chrom"], r["start"], r["end"], r["strand"],
                r["role"], r["link_internal"], r["link_reason"],
                r["u3_end_rel"], r["u5_start_rel"], r["confidence"], r["status"], r["evidence"]
            ]))+"\n")

    with bed_ok_path.open("w") as f:
        for (c,s,e,n,sc,st) in segs_ok:
            f.write(f"{c}\t{s}\t{e}\t{n}\t{sc}\t{st}\n")

    with bed_all_path.open("w") as f:
        for (c,s,e,n,sc,st) in segs_all:
            f.write(f"{c}\t{s}\t{e}\t{n}\t{sc}\t{st}\n")

    with flank_path.open("w") as f:
        for (c,s,e,n,sc,st) in flank:
            f.write(f"{c}\t{s}\t{e}\t{n}\t{sc}\t{st}\n")

    with signals_path.open("w") as f:
        for (c,s,e,n,sc,st) in signals:
            f.write(f"{c}\t{s}\t{e}\t{n}\t{sc}\t{st}\n")

    # Summary
    status_counts = Counter(r["status"] for r in rows)
    class_counts  = Counter(class_for_name(r["name"]) or "NA" for r in rows)
    role_counts = Counter(name2role.values())
    print(f"[OK] Components map:     {map_path}", file=sys.stderr)
    print(f"[OK] Components BED:     {bed_components_path}", file=sys.stderr)
    print(f"[OK] Catalogue:          {tsv_path}", file=sys.stderr)
    print(f"[OK] Segments (OK):      {bed_ok_path}", file=sys.stderr)
    print(f"[OK] Segments (all):     {bed_all_path}", file=sys.stderr)
    print(f"[OK] Flanks:             {flank_path}", file=sys.stderr)
    print(f"[OK] Signals:            {signals_path}", file=sys.stderr)
    print(f"[LOG] {logf}", file=sys.stderr)

    print("Summary (roles):", file=log)
    for k,v in role_counts.most_common(): print(f"  {k:10s} {v}", file=log)
    print("Summary (status):", file=log)
    for k,v in status_counts.most_common(): print(f"  {k:10s} {v}", file=log)
    print("Summary (classes):", file=log)
    for k,v in class_counts.most_common(): print(f"  {k:10s} {v}", file=log)
    log.close()




if __name__ == "__main__":
    main()
