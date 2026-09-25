#!/bin/bash
#SBATCH --job-name=bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100G

cd ./PREDICTION
# Run the Nextflow main pipeline
nextflow run .PREDICTION/main.nf \
    --input_csv ./path/to/input.csv \
    --outdir ./path/to/output \
    --reference_dir /path/to/reference \
    --gender_threshold \
    --workdir ./work \
    -with-trace ./path/to/output/trace_prediction.tsv \
    -with-timeline ./path/to/output/timeline_prediction.html \
    -with-report ./path/to/output/report_prediction.html \
    2>&1 | tee "./path/to/output/prediction.log"
