# BBETD: Blast-Based Extended Tandem duplication gene Detector

[![Language: Python](https://img.shields.io/badge/Language-Python%203.8+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**BBETD** is an efficient tool designed to identify tandemly duplicated genes (TDGs) across eukaryotic genomes. It leverages protein sequence similarity (Diamond with BLOSUM45) combined with gene synteny and a dual-rule classification system to identify both immediate (k=0) and extended (intervening genes, k ≤ 5) tandem arrays.

---

##  Installation

### Using Conda (Recommended)

```bash
git clone https://github.com/<your_username>/BBETD.git
cd BBETD

# Create environment and install dependencies (including Diamond & gffread)
conda env create -f environment.yml
conda activate bbetd

# Install BBETD command line
pip install -e .
```

---

##  Quick Start

```bash
bbetd \
  -g genome.fa \
  -a annotation.gff3 \
  -o ./results \
  -p my_species \
  -t 8
```

### Outputs

The output directory contains:
- `*.pTpG.pep`: Cleaned longest-isoform protein FASTA.
- `*.Order.tsv`: Gene synteny ordering table.
- `*_tandem_blocks.tsv`: Clustered tandem array blocks.
- `*_bbetd_tdg.pairs`: **Final high-confidence pairwise tandem duplicate gene list.**

---

##  Parameters

| Option | Default | Description |
|---|---|---|
| `-g, --genome` | *Required* | Path to reference genome FASTA |
| `-a, --gff` | *Required* | Path to annotation GFF3 / GTF file |
| `-o, --outdir` | `./bbetd_results` | Output directory path |
| `-p, --prefix` | `bbetd` | Prefix for output files |
| `-t, --threads` | `8` | Number of CPU threads |
| `--max-gap` | `5` | Maximum intervening genes allowed ($k \le 5$) |
| `--evalue` | `1e-2` | Diamond BLAST E-value cut-off |
| `--min-cov` | `20.0` | Minimum sequence coverage (%) |
| `--matrix` | `BLOSUM45` | Diamond alignment matrix |


##  Demo

```bash
cd scripts
bash run_slurm_example.sh
#[1/4] Extracting proteins & calculating gene orders...
#[2/4] Running Diamond self-blastp alignment..
# ...
#Reported 59'778 pairwise alignments, 59'778 HSPs.
#8'047 queries aligned.
#[3/4] Detecting tandem duplication blocks...
#[4/4] Performing Dual-Rule classification...
#==========================================================================
#[✓] Completed in 20.84s!
#[✓] High-confidence tandem duplicate pairs found: 1252
#[✓] Final results saved to: ./results/simulated_bbetd_bbetd_tdg.pairs
#==========================================================================
```
