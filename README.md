

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

## 📋 Example Usage

### Step 1: Extract LTR coordinates from RepeatMasker

This step parses the RepeatMasker .out file and outputs a BED file containing only LTR/ERV elements (excluding internal regions).

```
python 01_extract_ltr_coordinates_from_repeatmasker.py \
  --repeatmasker GRCh38.primary_assembly.genome.fa.out \
  --output results/LTR_raw.bed
```

### Step 2: Merge fragmented LTR annotations

RepeatMasker often splits LTRs into multiple fragments. This step merges nearby fragments by subfamily and strand using a 100 bp window by default. You can modify the `ALLOWED_GAP` parameter inside the script.

``` 
python 02_merge_close_ltrs_by_subfamily.py \
  results/LTR_raw.bed \
  results/LTR_merged.bed
```

### Step 3: Extract LTR genomic sequences

Using the merged BED file and the reference genome, this step produces a strand-corrected FASTA file of LTRs.

```
python 03_extract_ltr_sequences.py \
  --bed results/LTR_merged.bed \
  --genome GRCh38.primary_assembly.genome.fa \
  --output results/LTR_merged.fa
```

### Step 4: Split LTR FASTA by subfamily

Each subfamily is written to its own FASTA file, which will be used for FIMO scanning.

```
python 04_split_fasta_by_subfamily.py \
  --input results/LTR_merged.fa \
  --output_dir results/by_subfamily
```

### Step 5: Build a global 0th-order background model

``` 
bash 05_fimo_create_background.sh \
  results/LTR_merged.fa \
  results/bg_LTR_all.txt
```

### Step 6: Run FIMO motif scanning

This step scans each subfamily FASTA using JASPAR 2024 motifs (converted to MEME format).

``` 
for fa in results/by_subfamily/*.fa; do
    sub=$(basename "$fa" .fa)
    python 06_run_fimo.py \
      --motifs JASPAR2024_CORE_vertebrates.meme \
      --fasta "$fa" \
      --bgfile results/bg_LTR_all.txt \
      --output results/fimo/$sub
done
``` 

### Step 7: Merge all FIMO outputs

Combine all per-subfamily fimo.tsv tables into one file.

```
python 07_merge_all_fimo.py \
  --input-dir results/fimo \
  --output results/merged_fimo.tsv
```

### Step 8: Parse FIMO hits and generate BED coordinates

This step converts motif hits into genomic coordinates and outputs two files:

* sorted TSV with all hits
* BED file containing all motifs coordinates

```
python 08_parse_fimo_and_generate_bed.py \
  --input results/merged_fimo.tsv \
  --output_tsv results/merged_fimo_sorted.tsv \
  --output_bed results/merged_fimo.bed
```

### Step 9: Build a PBS tRNA mini-library

This step extracts 3′ tRNA tails from GtRNAdb and prepares DNA and reverse-complement libraries for PBS detection.

```
python 09_build_trna_minilib.py \
  --input hg38-mature-tRNAs.fa \
  --outdir motifs/tRNA_PBS \
  --tail 20
``` 

### Step 10: Annotate U3–R–U5 segments and flanks (PBS/PPT)

This is the core step that assigns each LTR as 5′/3′/solo, detects promoter elements, PAS, PBS, and PPT, and outputs segmentation BED/TSVs.

```
python 10_annotate_u3r_u5.py \
  --ltr-bed results/LTR_merged.bed \
  --ltr-fasta results/LTR_merged.fa \
  --genome GRCh38.primary_assembly.genome.fa \
  --internal-bed results/internal_regions.bed \
  --subfamily-map repeats_classfamily_map.tsv \
  --promoters promoter_motifs.tsv \
  --pas pas_hexamers.txt \
  --trna-rc motifs/tRNA_PBS/tRNA_3prime_DNA_revcomp_dedup.fa \
  --outdir results/U3R_U5
``` 

### Step 11: Downstream analysis (optional)

To reproduce the plots and analyses used in the manuscript, you can run the accompanying R scripts:

``` 
Rscript analysis_ltr_regulation.R
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

🔗 [https://doi.org/10.5281/zenodo.17602210](https://doi.org/10.5281/zenodo.17602210)

Please cite our article (see below) and the corresponding dataset DOI if you reuse the data in your own work.

## 🧠 Citation

If you use this resource, please cite:

> Montserrat-Ayuso, T., & Esteve-Codina, A. (2025). *Regulatory Features and Functional Specialization of Human Endogenous Retroviral LTRs: A Genome-Wide Annotation and Analysis via HERVarium.*
> *bioRxiv* (in preparation).

## 📎 License

MIT License. See LICENSE for details.
