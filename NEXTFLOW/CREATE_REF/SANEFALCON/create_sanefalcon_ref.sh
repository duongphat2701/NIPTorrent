cd /mnt/d18t/DANPHAM/NIPT/NIPTorrent/CREATE_REF
# Run the Nextflow pipeline
## SIZE 454
# nextflow run /mnt/d18t/DANPHAM/NIPT/NIPTorrent/CREATE_REF/main1.nf \
#     --input_csv /mnt/d18t/DANPHAM/NIPT/NIPT_DATA/REF/TRIMMING_FASTQ/TRIM_15_50/mapped/sample.csv \
#     --outdir /mnt/d18t/DANPHAM/NIPT/NIPTorrent/CREATE_REF/SANEFALCON/REF \
#     --ff_ref_tsv /mnt/d18t/DANPHAM/NIPT/NIPTorrent/CREATE_REF/SANEFALCON/ref_ff_454.tsv \
#     -work-dir /mnt/d18t/DANPHAM/NIPT/NIPTorrent/CREATE_REF/SANEFALCON/work 

export SANEFALCON_DIR="/mnt/d18t/DANPHAM/NIPT/NIPTorrent/CREATE_REF/assets/sanefalcon"
export FF_REF_TSV="/mnt/d18t/DANPHAM/NIPT/NIPTorrent/CREATE_REF/SANEFALCON/ff_ref.tsv"
export OUTDIR="/mnt/d18t/DANPHAM/NIPT/NIPTorrent/CREATE_REF/SANEFALCON/REF"

bash /mnt/d18t/DANPHAM/NIPT/NIPTorrent/CREATE_REF/bin/sanefalcon.sh \
    --sanefalcon-dir $SANEFALCON_DIR \
    --bamdir "/mnt/d18t/DANPHAM/NIPT/NIPT_DATA/REF/TRIMMING_FASTQ/TRIM_15_50/mapped" \
    --outdir $OUTDIR \
    --threads 32 \
    --ff_ref_tsv $FF_REF_TSV \
    --all