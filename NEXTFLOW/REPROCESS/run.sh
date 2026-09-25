#!/bin/bash
#SBATCH --job-name=NIPT_REPROCESS
#SBATCH --output=NIPT_REPROCESS_%j.out
#SBATCH --error=NIPT_REPROCESS_%j.err
#SBATCH --cpus-per-task=40
#SBATCH --mem=100G

cd ./NIPTorrent/REPROCESS
# Run the Nextflow main pipeline
#    --trimmomatic_options MINLEN:50 
#    --trimmomatic_options MINLEN:50 SLIDINGWINDOW:30:15
nextflow ./NIPTorrent/REPROCESS/main.nf \
    --input_csv ./path/to/input.csv \
    --outdir ./path/to/output \
    -with-trace ./path/to/output/trace1.tsv \
    -with-timeline ./path/to/output/timeline1.html \
    -with-report ./path/to/output/report1.html \
    --workDir ./path/to/output/work
