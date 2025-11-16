

# HERV Regulatory Annotation

## Genome-wide U3–R–U5 segmentation and regulatory motif mapping of HERV LTRs
This repository contains the complete pipeline for reconstructing and annotating the regulatory architecture of human endogenous retrovirus (HERV) long terminal repeats (LTRs). It includes scripts for merging fragmented RepeatMasker elements, assigning LTRs to internal ERV regions, delineating U3–R–U5 segments in a strand-aware manner, scanning regulatory motifs (TFBMs), and identifying flanking PBS/PPT sequences.  
<p align="center">
 <img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/247104a8-3b0f-4023-926d-7032fb187212" />
</p>

## 🧬 Overview
LTRs encode the *regulatory interface* of retroviruses. Although most human HERVs are ancient and degenerated, many retain recognizable promoters, transcription-factor binding motifs (TFBMs), polyadenylation signals, or tRNA/PPT features that once drove retroviral replication and now may shape host gene regulation.

Here, we provide a unified, genome-wide framework for:

* reconstructing HERV LTR structures
* segmenting each LTR into U3, R, and U5 regions
* classifying LTRs as 5′, 3′, tandem, or solo
* detecting PBS and PPT flanking signals
* scanning >750 vertebrate TF motifs using FIMO
* quantifying TFBM burden and positional enrichment
* generating BED/TSV files for downstream analyses

This project forms the second major component of **HERVarium**, alongside the internal domain annotation (*herv-domain-map*).

## 🚀 Features

* RepeatMasker-based LTR reconstruction:

  * Extract LTR coordinates directly from RepeatMasker .out files
  * Merge adjacent or fragmented LTR pieces by subfamily and strand, using a distance threshold
  * Recover LTR genomic sequences from the reference genome FASTA
  * Split merged LTRs into per-subfamily FASTA files for downstream FIMO scanning

* Motif scanning (TFBM detection)

  * Build a single, global 0th-order Markov background from all LTR sequences
  * Run FIMO (MEME suite) for each subfamily using JASPAR 2024 PFMs converted to MEME format
  * Merge all FIMO results across subfamilies into a standardized table

* U3–R–U5 structural segmentation

  * Assign each LTR as 5′, 3′ or solo based on proximity to internal proviral regions
  * Detect promoter-associated motifs within U3 (e.g., TATA box, Inr, BRE, DPE)
  * Identify PAS (polyadenylation signal) hexamers in U3/R/U5
  * Detect canonical PBS (tRNA-binding site) downstream of 5′ LTRs and PPT (polypurine tract) upstream of 3′ LTRs
  * Output U3–R–U5 segments with confidence scores and problem-flagging

## 🛠️ Dependencies

Conda environment:
```
name: ervreg
channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - python=3.8
  - pip
  - pyfaidx=0.8.1.4
  - biopython=1.83
  - pandas=2.0.3
  - pysam=0.22.1
  - tqdm
  - meme
  - samtools
  - bedtools
```


## 📁 Repository Structure
```
.
├── scripts/
│   ├── 01_extract_ltr_coordinates_from_repeatmasker.py
│   ├── 02_merge_close_ltrs_by_subfamily.py
│   ├── 03_extract_ltr_sequences.py
│   ├── 04_split_fasta_by_subfamily.py
│   ├── 05_fimo_create_background.sh
│   ├── 06_run_fimo.py
│   ├── 07_merge_all_fimo.py
│   ├── 08_parse_fimo_and_generate_bed.py
│   ├── 09_generate_ltr_catalogue.py
│   ├── 09_build_trna_minilib.py
│   ├── 10_annotate_u3r_u5.py
├── README.md
├── LICENSE
```

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

🔗 [https://doi.org/10.5281/zenodo.16318928](https://doi.org/10.5281/zenodo.16318928)

Please cite our article (see below) and the corresponding dataset DOI if you reuse the data in your own work.

## 🧠 Citation

If you use this resource, please cite:

> Montserrat-Ayuso, T., & Esteve-Codina, A. (2025). *Regulatory Features and Functional Specialization of Human Endogenous Retroviral LTRs: A Genome-Wide Annotation and Analysis via HERVarium.*
> *bioRxiv* (in preparation).

## 📎 License

MIT License. See LICENSE for details.
