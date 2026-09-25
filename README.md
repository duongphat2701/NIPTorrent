# NIPTorrent: NIPT Data Analysis for the Ion Torrent Sequencing Platform

Non-invasive prenatal testing (NIPT) is a crucial screening tool for fetal chromosomal aneuploidies, leveraging cell-free fetal DNA (cffDNA) in maternal blood. While the Ion Torrent sequencing platform—with its cost-effective semiconductor technology and Torrent Suite software—presents a promising alternative for NIPT, there is a significant gap: the absence of a published, open-source bioinformatics pipeline specifically tailored for it. 

This gap hinders accurate sequencing data processing, fetal fraction (ff) estimation, and z-score calculation, which are essential for reliable aneuploidy detection. 

To address this, we developed **NIPTorrent**, a comprehensive bioinformatics pipeline built with Nextflow. We constructed a Vietnamese-specific reference dataset and benchmarked key computational methods for ff estimation and z-score calculation to optimize their performance on Ion Torrent data. Our pipeline enhances accuracy, sensitivity, and data processing efficiency by refining fetal fraction assessment, statistical modeling, and aneuploidy detection. This validated, reproducible solution facilitates the clinical integration of Ion Torrent sequencing, contributing to expanded access to cost-effective prenatal screening and improved early fetal risk assessment.

**Keywords:** Non-Invasive Prenatal Testing (NIPT), cell-free fetal DNA (cffDNA), bioinformatics, Torrent Suite, Ion Torrent sequencing platform and open-source bioinformatics pipeline

---

## 1. Input Data Preparation

The pipeline requires a standard CSV file to map your single-ended sequencing data. 

### CSV Format Requirements

| Column | Description |
| :--- | :--- |
| `sample` | Unique identifier for the sample |
| `bam` | Absolute path to the sample's `.bam` file |
| `bai` | Absolute path to the sample's `.bai` index file |

**Example `input.csv`:**
```csv
sample,bam,bai
sample1,/path/to/sample1.bam,/path/to/sample1.bai
sample2,/path/to/sample2.bam,/path/to/sample2.bai
```

### Automated CSV Generation

You can use the following helper script to automatically generate the required CSV file from a directory containing your BAM files.

```bash
MAP_dir="path/to/sample_dir"
output="path/to/file.csv"

# Create header
echo "sample,bam,bai" > "$output"

# Loop through all BAM files and append to CSV
for bam in "$MAP_dir"/*.bam; do
    sample_name=$(basename "$bam" .bam)
    bai="$MAP_dir/${sample_name}.bam.bai"
    echo "$sample_name,$bam,$bai" >> "$output"
done
```

---

## 2. Reference Datasets & Thresholds

We provide a reference set built from **500 normal samples**. Depending on the read trimming condition applied, the total number of remaining unique reads changes, which consequently alters the required gender threshold. 

Choose the appropriate reference path and corresponding gender threshold for your analysis from the table below:

| Condition | Trimming Parameters | Gender Threshold | Reference Directory Path |
| :--- | :--- | :--- | :--- |
| **Raw** | Non-trimming (Raw) | 0.03 | `/mnt/d18t/DANPHAM/NIPT/NIPTorrent/PREDICTION/REFERENCES/REF_500_old/RAW` |
| **Length_only** | Length: 50 bp | 0.03 | `/mnt/d18t/DANPHAM/NIPT/NIPTorrent/PREDICTION/REFERENCES/REF_500_old/TRIM_50` |
| **Length_quality** | Length: 50 bp, Quality: 15 | 0.02 | `/mnt/d18t/DANPHAM/NIPT/NIPTorrent/PREDICTION/REFERENCES/REF_500_old/TRIM_15_50` |
| **Length_quality** | Length: 50 bp, Quality: 20 | 0.02 | `/mnt/d18t/DANPHAM/NIPT/NIPTorrent/PREDICTION/REFERENCES/REF_500_old/TRIM_20_50` |

> **Note:** Please ensure you fill in the missing gender threshold for the Quality 20 trimming condition before finalizing this documentation.

---

## 3. Running the Pipeline

Execute the main Nextflow pipeline using the following command. Ensure you replace the placeholder paths and threshold values with the parameters specific to your current run.

```bash
cd /mnt/d18t/DANPHAM/NIPT/NIPTorrent/PREDICTION

nextflow run main.nf \
    --input_csv /path/to/csv_files.csv \
    --outdir /path/to/destination_directory \
    --reference_dir /path/to/references_directory \
    --gender_threshold [insert_number_from_table] \
    --workdir /mnt/d18t/DANPHAM/NIPT/NIPTorrent/PREDICTION/work \
    -with-trace /path/to/destination_directory/trace_prediction.tsv \
    -with-timeline /path/to/destination_directory/timeline_prediction.html \
    -with-report /path/to/destination_directory/report_prediction.html \
    2>&1 | tee "/path/to/destination_directory/prediction.log"
```