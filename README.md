# NIPTorrent: NIPT Data Analysis for the Ion Torrent Sequencing Platform

**Non-Invasive Prenatal Testing (NIPT) bioinformatics pipeline for Ion Torrent sequencing data.**

---

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Pipeline 1: REPROCESS](#pipeline-1-reprocess)
4. [Pipeline 2: CREATE_REF](#pipeline-2-create_ref)
5. [Pipeline 3: PREDICTION](#pipeline-3-prediction)
6. [Input File Format](#input-file-format)
7. [Reference Datasets](#reference-datasets)
8. [Quick Start](#quick-start)

---

## Overview

**NIPTorrent** is a comprehensive bioinformatics pipeline for analyzing Non-Invasive Prenatal Testing (NIPT) data generated on the Ion Torrent sequencing platform. The pipeline processes raw sequencing data through quality control, alignment, and advanced statistical analyses to detect fetal chromosomal aneuploidies.

### Key Features

- **Multi-mapper support**: BWA-MEM, Bowtie2, and TMAP aligners
- **Fetal Fraction Estimation**: Three methods (SeqFF, DeFrag, SaneFALCON)
- **Aneuploidy Detection**: WisecondorX and NIPTeR algorithms
- **Reproducible**: Built with Nextflow for containerized, scalable execution
- **Optimized for Ion Torrent**: Benchmarked specifically for this sequencing platform

---

## Pipeline Architecture

The NIPTorrent workflow consists of **3 independent pipelines** that must be executed in sequence:

```
┌─────────────────────────────────────────────────────────────────┐
│                        INPUT DATA                                │
│                   (Unmapped BAM files)                          │
└─────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│              PIPELINE 1: REPROCESS                               │
│  UBAM → FASTQ → Mapping → Sort → Remove Duplicates → Index    │
│  (BWA-MEM / Bowtie2 / TMAP)                                    │
└─────────────────────────────────────────────────────────────────┘
                    │                       │
                    ▼                       ▼
        ┌───────────────────┐   ┌───────────────────┐
        │   PIPELINE 2:     │   │   PIPELINE 3:     │
        │   CREATE_REF       │   │   PREDICTION      │
        │                   │   │                   │
        │   Build reference │   │   Fetal fraction  │
        │   datasets from   │   │   Aneuploidy      │
        │   normal samples  │   │   detection       │
        └───────────────────┘   └───────────────────┘
```

### Step-by-step Process Flow

#### REPROCESS Pipeline

| Step | Tool | Description |
|------|------|-------------|
| 1 | **samtools fastq** | Convert unmapped BAM (UBAM) to FASTQ format |
| 2 | **Trimmomatic** | Quality trimming (optional, can be skipped) |
| 3 | **FastQC** | Pre/Post quality control reports |
| 4 | **MultiQC** | Aggregate QC reports |
| 5 | **BWA-MEM** | Map reads to reference genome |
| 6 | **Bowtie2** | Alternative mapper (very-sensitive mode) |
| 7 | **TMAP** | Ion Torrent-specific mapper |
| 8 | **samtools sort** | Sort BAM files by coordinate |
| 9 | **samtools markdup** | Remove duplicate reads |
| 10 | **samtools view -F 4** | Filter to keep only mapped reads |
| 11 | **samtools index** | Index final BAM files |

#### CREATE_REF Pipeline

Creates reference datasets required for PREDICTION pipeline:

| Tool | Description |
|------|-------------|
| **SaneFALCON** | Build fetal fraction reference |
| **NIPTeR** | Create control group reference |
| **convert_files** | Generate .gcc, .pickle, .npz files |
| **convert_gender** | Create gender-specific references |

#### PREDICTION Pipeline

| Analysis | Tool | Description |
|----------|------|-------------|
| **Fetal Fraction** | SeqFF | Fast estimation from read counts |
| | DeFrag | Fragment-based estimation (2 variants) |
| | SaneFALCON | GC-corrected estimation |
| **Aneuploidy Detection** | WisecondorX | Z-score based detection with/without blacklist masking |
| | NIPTeR | NCV (Normalized Chromosome Value) based detection |
| **Quality Control** | PCA Plot | Sample quality visualization |
| | Gender Prediction | Determine fetal/maternal sample gender |

---

## Pipeline 1: REPROCESS

Converts unmapped BAM files to analysis-ready sorted, deduplicated BAM files using multiple aligners.

### Directory Structure

```
NEXTFLOW/REPROCESS/
├── main.nf                    # Main pipeline
├── run.sh                     # Execution script
├── UPSTREAM/
│   ├── bwa.nf                 # BWA-MEM alignment
│   ├── bowtie2.nf             # Bowtie2 alignment
│   ├── tmap.nf                # TMAP alignment (Ion Torrent)
│   └── samtools/
│       ├── sort_bam.nf        # Sort BAM
│       ├── rm_dup.nf          # Remove duplicates
│       ├── mapped_bam.nf       # Filter mapped reads
│       └── index.nf           # Index BAM
```

### Parameters

| Parameter | Description | Required |
|-----------|-------------|----------|
| `--input_csv` | Path to input CSV file | Yes |
| `--outdir` | Output directory | Yes |
| `--workDir` | Nextflow work directory | Yes |
| `--ref_bwa` | BWA reference genome index | Yes |
| `--ref_bowtie2` | Bowtie2 reference index | Yes |
| `--ref_tmap` | TMAP reference | Yes |
| `--trimmomatic_options` | Trimmomatic parameters | No |

### Example Command

```bash
cd NEXTFLOW/REPROCESS

nextflow run main.nf \
    --input_csv /path/to/samples.csv \
    --outdir /path/to/output \
    --ref_bwa /path/to/hg38.fa \
    --ref_bowtie2 /path/to/hg38 \
    --ref_tmap /path/to/tmap_ref \
    --workDir /path/to/output/work \
    -with-trace /path/to/output/trace.tsv \
    -with-timeline /path/to/output/timeline.html \
    -with-report /path/to/output/report.html
```

### Output Directory Structure

```
outdir/
├── sorted_bwa/          # BWA-MEM sorted BAM files
│   ├── sample1_sorted.bam
│   └── sample1_sorted.bam.bai
├── dedup_bwa/           # Deduplicated BAM files
├── mapped_bwa/          # Mapped-only BAM files (unmapped removed)
├── sorted_bowtie2/      # Bowtie2 results
├── dedup_bowtie2/
├── mapped_bowtie2/
├── sorted_tmap/         # TMAP results
├── dedup_tmap/
└── mapped_tmap/
```

---

## Pipeline 2: CREATE_REF

Creates reference datasets from normal (euploid) samples for use in the PREDICTION pipeline.

### Directory Structure

```
NEXTFLOW/CREATE_REF/
├── main.nf                    # Main pipeline
├── run.sh                     # Execution script
├── DOWNSTREAM/
│   ├── create_ref_sanefalcon.nf   # SaneFALCON reference
│   ├── create_ref_nipter.nf       # NIPTeR reference
│   ├── create_ref_nipter1.nf      # NIPTeR reference (new)
│   ├── convert_files.nf           # File format conversion
│   └── convert_gender.nf         # Gender reference
├── bin/
│   └── sanefalcon.sh              # SaneFALCON wrapper
└── assets/sanefalcon/            # SaneFALCON scripts
```

### Parameters

| Parameter | Description | Required |
|-----------|-------------|----------|
| `--input_csv` | Path to input CSV file with normal samples | Yes |
| `--binSize` | Bin size for genome binning (default: 1000000) | No |
| `--outdir` | Output directory | Yes |
| `--workDir` | Nextflow work directory | Yes |
| `--ff_ref_tsv` | Fetal fraction reference TSV file | Yes (for SaneFALCON) |

### Example Command

```bash
cd NEXTFLOW/CREATE_REF

nextflow run main.nf \
    --input_csv /path/to/normal_samples.csv \
    --binSize 1000000 \
    --outdir /path/to/reference_output \
    --ff_ref_tsv /path/to/ff_reference.tsv \
    --workDir /path/to/reference_output/work \
    -with-trace /path/to/reference_output/trace.tsv \
    -with-timeline /path/to/reference_output/timeline.html \
    -with-report /path/to/reference_output/report.html
```

### Output Directory Structure

```
outdir/
├── sanefalcon_ref/            # SaneFALCON reference
│   ├── profile_ref.tsv
│   └── model files...
├── NIPTeR/reference/          # NIPTeR control group
│   ├── NIPTeR_ref.controlgroup.GC.rds
│   └── NIPTeR_ref.controlgroup.diagnostics.rds
├── converted_files/           # Converted sample files
│   ├── sample1.gcc
│   ├── sample1.pickle
│   └── sample1.npz
└── gender/                    # Gender prediction references
```

> **Important**: The output from CREATE_REF should be used as `--reference_dir` in the PREDICTION pipeline.

---

## Pipeline 3: PREDICTION

Analyzes samples for fetal fraction estimation and chromosomal aneuploidy detection.

### Directory Structure

```
NEXTFLOW/PREDICTION/
├── main.nf                        # Main pipeline
├── run.sh                         # Execution script
├── DOWNSTREAM/
│   ├── seqff.nf                   # SeqFF fetal fraction
│   ├── defrag_a.nf                # DeFrag method A
│   ├── defrag_b.nf                # DeFrag method B
│   ├── wisecondorX.nf             # WisecondorX detection
│   ├── nipter.nf                  # NIPTeR detection
│   ├── blacklist_size_mask.nf     # Blacklist generation
│   ├── predict_abnormal_*.nf       # Prediction with/without blacklist
│   ├── convert_files.nf            # File format conversion
│   ├── convert_gender.nf          # Gender prediction
│   ├── read_count.nf              # Read counting
│   ├── PCA_plot.nf                # PCA visualization
│   ├── combine.nf                  # Combine results
│   └── report_NIPT.nf              # Final report
└── bin/
    ├── seqff.py
    ├── wisecondorX.py
    ├── nipter1.R
    └── ...
```

### Parameters

| Parameter | Description | Required |
|-----------|-------------|----------|
| `--input_csv` | Path to input CSV file | Yes |
| `--outdir` | Output directory | Yes |
| `--reference_dir` | Directory with reference datasets | Yes |
| `--gender_threshold` | Threshold for gender prediction | Yes |
| `--basic_blacklist` | Path to blacklist BED file | No |
| `--binSize` | Bin size (default: 1000000) | No |
| `--workDir` | Nextflow work directory | Yes |

### Example Command

```bash
cd NEXTFLOW/PREDICTION

nextflow run main.nf \
    --input_csv /path/to/samples.csv \
    --outdir /path/to/results \
    --reference_dir /path/to/reference_output \
    --gender_threshold 0.03 \
    --workDir /path/to/results/work \
    -with-trace /path/to/results/trace.tsv \
    -with-timeline /path/to/results/timeline.html \
    -with-report /path/to/results/report.html \
    2>&1 | tee /path/to/results/prediction.log
```

### Output Directory Structure

```
outdir/
├── converted_files/               # Converted BAM files
├── seqff/                         # SeqFF results
│   └── seqff_summary.csv
├── defrag_a/                      # DeFrag A results
├── defrag_b/                      # DeFrag B results
│   └── defrag_b_summary.csv
├── wisecondorX/                   # WisecondorX results
│   ├── aberrations/
│   ├── z_bins/
│   ├── z_segments/
│   ├── statistics/
│   ├── abnormal_tables/
│   └── wisecondorX_summary.csv
├── abnormal_without_blacklist_prediction/
├── abnormal_size_mask_blacklist/
└── NIPTeR/                        # NIPTeR results (if enabled)
    ├── rds/
    └── abnormal_tables/
```

---

## Input File Format

All pipelines require a CSV file with the following format:

### CSV Format

```csv
sample,bam,bai
sample1,/path/to/sample1.bam,/path/to/sample1.bai
sample2,/path/to/sample2.bam,/path/to/sample2.bai
```

| Column | Description |
|--------|-------------|
| `sample` | Unique sample identifier |
| `bam` | Absolute path to BAM file |
| `bai` | Absolute path to BAI index file |

### Helper Script to Generate CSV

```bash
#!/bin/bash
MAP_dir="path/to/bam_directory"
output="path/to/samples.csv"

# Create header
echo "sample,bam,bai" > "$output"

# Loop through all BAM files and append to CSV
for bam in "$MAP_dir"/*.bam; do
    sample_name=$(basename "$bam" .bam)
    bai="$bam.bai"
    echo "$sample_name,$bam,$bai" >> "$output"
done

echo "CSV file created: $output"
```

---

## Reference Datasets

Pre-built reference datasets are available for Vietnamese population. These were generated from 500 normal samples.

### Available References

| Condition | Trimming Parameters | Gender Threshold | Reference Directory |
|-----------|---------------------|-----------------|---------------------|
| **Raw** | None | 0.03 | `REF_500_old/RAW` |
| **Trim_50** | Length: 50 bp | 0.03 | `REF_500_old/TRIM_50` |
| **Trim_15_50** | Length: 50 bp, Quality: 15 | 0.02 | `REF_500_old/TRIM_15_50` |
| **Trim_20_50** | Length: 50 bp, Quality: 20 | 0.02 | `REF_500_old/TRIM_20_50` |

### Reference Directory Structure

```
reference_dir/
├── boydir/                    # Male-specific bins
├── girldir/                   # Female-specific bins
├── wisecondorx_reference.npz # WisecondorX reference
├── NIPTeR_control.rds        # NIPTeR control group
├── sanefalcon/               # SaneFALCON profile
│   ├── profile_ref.tsv
│   └── ...
└── ff_reference.tsv          # Fetal fraction reference
```

---

## Quick Start

### Complete Workflow Example

```bash
# 1. Run REPROCESS pipeline
cd NEXTFLOW/REPROCESS
nextflow run main.nf \
    --input_csv /data/samples.csv \
    --outdir /results/reprocessed \
    --ref_bwa /ref/hg38.fa \
    --ref_bowtie2 /ref/hg38 \
    --ref_tmap /ref/tmap_ref \
    --workDir /results/reprocessed/work

# 2. Create reference from normal samples (first time only)
cd ../CREATE_REF
nextflow run main.nf \
    --input_csv /data/normal_samples.csv \
    --binSize 1000000 \
    --outdir /reference/my_reference \
    --ff_ref_tsv /reference/ff_ref.tsv \
    --workDir /reference/my_reference/work

# 3. Run PREDICTION pipeline
cd ../PREDICTION
nextflow run main.nf \
    --input_csv /results/reprocessed/mapped_bwa/samples.csv \
    --outdir /results/prediction \
    --reference_dir /reference/my_reference \
    --gender_threshold 0.03 \
    --workDir /results/prediction/work
```

---

## Support

For issues or questions, please refer to the project documentation or contact the maintainers.

---

**Keywords:** Non-Invasive Prenatal Testing (NIPT), cell-free fetal DNA (cffDNA), bioinformatics, Ion Torrent, aneuploidy detection, fetal fraction estimation
