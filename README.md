

# HERV Regulatory Annotation

## Genome-wide U3–R–U5 segmentation, regulatory motif mapping, and structural annotation of HERV LTRs
This repository contains the complete pipeline for reconstructing and annotating the regulatory architecture of human endogenous retrovirus (HERV) long terminal repeats (LTRs). It includes scripts for merging fragmented RepeatMasker elements, assigning LTRs to internal ERV regions, delineating U3–R–U5 segments in a strand-aware manner, scanning regulatory motifs (TFBMs), and identifying flanking PBS/PPT sequences.

---

## 🧬 Overview
LTRs encode the *regulatory interface* of retroviruses. Although most human HERVs are ancient and degenerated, many retain recognizable promoters, transcription-factor binding motifs (TFBMs), polyadenylation signals, or tRNA/PPT features that once drove retroviral replication and now may shape host gene regulation.

herv-regulatory-map provides a unified, genome-wide framework for:

* reconstructing HERV LTR structures
* segmenting each LTR into U3, R, and U5 regions
* classifying LTRs as 5′, 3′, tandem, or solo
* detecting PBS and PPT flanking signals
* scanning >750 vertebrate TF motifs using FIMO
* quantifying TFBM burden and positional enrichment
* generating publication-ready BED/TSV files for downstream analyses

This project forms the second major component of **HERVarium**, alongside the internal domain annotation (*herv-domain-map*).

---

## 🚀 Features

* RepeatMasker-based LTR reconstruction

  * Merge fragmented LTRs
  * Associate LTRs with nearby internal regions
  * Classify as 5′, 3′, both (tandem), or solo

* U3–R–U5 segmentation

  * Strand-aware reconstruction
  * Boundary discovery and confidence scoring
  * High- and low-confidence outputs

* Regulatory signal annotation

  * TFBM scanning using JASPAR 2024 CORE (PFMs → MEME format)
  * Motif burden quantification (per LTR, per family)
  * Positional density analysis (KDE, histogram, peak detection)
  * Promoter-associated motifs
  * PAS hexamer detection
  * PBS (tRNA binding) and PPT detection

* BED/TSV output

  * Segmented LTRs
  * High-confidence vs all-confidence sets
  * PBS/PPT flanks
  * Motif hits (FIMO)
  * Summary tables per LTR, per family, and per regulatory category

---

## 🛠️ Dependencies

### Conda environment

```yaml
name: herv-regmap
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.10
  - bedtools
  - bioconda::pyfaidx
  - bioconda::meme
  - pandas
  - numpy
  - biopython
  - tqdm
```

Additional tools:

* **FIMO** (MEME Suite 5.x)
* **RepeatMasker** `.out` file for GRCh38
* **Genome FASTA (GRCh38 from GENCODE)**
* **JASPAR 2024** PFMs (CORE vertebrates, MEME format)

---

## 📁 Repository Structure
```
.
├── scripts/
│   ├── 01_merge_repeatmasker_ltrs.py
│   ├── 02_assign_ltr_roles.py
│   ├── 03_segment_u3_r_u5.py
│   ├── 04_detect_pbs_ppt.py
│   ├── 05_prepare_fasta_for_fimo.py
│   ├── 06_run_fimo.py
│   ├── 07_parse_fimo_results.py
│   ├── 08_summarize_regulatory_features.py
│   ├── 09_generate_ltr_catalogue.py
│   ├── run_ltr_annotation_pipeline.py
│   ├── run_ltr_motif_pipeline.py
│   ├── run_full_regulatory_pipeline.py
├── README.md
├── LICENSE
```

---

## 📋 Example Usage

### Step 1: LTR Reconstruction & U3–R–U5 Segmentation

```bash
python run_ltr_annotation_pipeline.py \
  --repeatmasker GRCh38.genome.fa.out \
  --genome GRCh38.genome.fa \
  --output results/LTR_annotation
```

This generates:

* `HERV_LTR_U3_R_U5_segments.bed` (high-confidence)
* `HERV_LTR_U3_R_U5_segments_all.bed` (all segments)
* `HERV_LTR_U3_R_U5_flanks.bed` (PBS/PPT)

### Step 2: Motif Scanning with FIMO

```bash
python run_ltr_motif_pipeline.py \
  --fasta results/LTR_annotation/LTRs.fa \
  --motifs JASPAR2024_CORE_vertebrates_non-redundant.meme \
  --bgfile results/bg_LTR_all.txt \
  --threshold 1e-4 \
  --output results/fimo_hits
```

### Step 3 Summary & Integration

```bash
python run_full_regulatory_pipeline.py \
  --segments results/LTR_annotation/HERV_LTR_U3_R_U5_segments.bed \
  --fimo results/fimo_hits \
  --output results/LTR_regulatory_summary
```
---

## 📦 Dataset

This repository contains the scripts used to generate the LTR regulatory dataset accompanying our study.

The full dataset—including:

* U3–R–U5 segment annotations
* LTR role classification (5′ / 3′ / tandem / solo)
* PBS/PPT sites
* TF motif burden tables
* Genome-wide FIMO hits
* Per-family regulatory summaries

is publicly available on Zenodo:

🔗 **[https://doi.org/10.5281/zenodo.16318928](https://doi.org/10.5281/zenodo.16318928)**

---

## 🧠 Citation

If you use this resource, please cite:

> **Montserrat-Ayuso, T., & Esteve-Codina, A. (2025). Regulatory architecture and transcription-factor motif landscape of human endogenous retrovirus LTRs.**
> *bioRxiv* (in preparation).
> Dataset: [https://doi.org/10.5281/zenodo.16318928](https://doi.org/10.5281/zenodo.16318928)

---

## 📎 **License**

MIT License. See LICENSE for details.
