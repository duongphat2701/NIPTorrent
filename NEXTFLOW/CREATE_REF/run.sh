#!/bin/bash
#SBATCH --job-name=NIPT_CREATE_REF
#SBATCH --output=NIPT_CREATE_REF_%j.out
#SBATCH --error=NIPT_CREATE_REF_%j.err
#SBATCH --cpus-per-task=40
#SBATCH --mem=100G

cd ./CREATE_REF
# Run the Nextflow main pipeline
nextflow ./CREATE_REF/main.nf \
    --input_csv ./path/to/input.csv \
    --binSize 1000000 \
    --outdir ./path/to/output \
    -with-trace ./path/to/output/trace.tsv \
    -with-timeline ./path/to/output/timeline.html \
    -with-report ./path/to/output/report.html \
    --workDir ./path/to/output/work