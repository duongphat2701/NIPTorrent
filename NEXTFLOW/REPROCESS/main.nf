#!/usr/bin/env nextflow

// Include upstream modules
include { convert_fastq } from "./UPSTREAM/samtools/fastq.nf"
include { fastqc_pre } from "./UPSTREAM/fastqc_pre.nf"
include { multiqc_pre } from "./UPSTREAM/multiqc_pre.nf"
include { trimmomatic } from "./UPSTREAM/trimmomatic.nf"
include { fastqc_post } from "./UPSTREAM/fastqc_post.nf"
include { multiqc_post } from "./UPSTREAM/multiqc_post.nf"
include { bwa } from "./UPSTREAM/bwa.nf"
include { sort_bam as sort_bwa } from "./UPSTREAM/samtools/sort_bam.nf"
include { rm_dup as rm_dup_bwa } from "./UPSTREAM/samtools/rm_dup.nf"
include { mapped_bam as mapped_bwa } from "./UPSTREAM/samtools/mapped_bam.nf"
include { index as index_bwa } from "./UPSTREAM/samtools/index.nf"
workflow {
    // Define input channels
    Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.read) }
        .set { input_channel }
    //--------------------------------------------
    // CONVERSION and TRIMMING
    //--------------------------------------------
    // Convert BAM to fastq.gz
    fastq = convert_fastq(input_channel)
    // Pre-trimming FastQC
    fastqc_pre = fastqc_pre(fastq)
    // MultiQC report before trimming
    multiqc_input_pre = fastqc_pre.collect()
    multiqc_pre = multiqc_pre(multiqc_input_pre)
    // Trimming with Trimmomatic
    trimmed_fastq = trimmomatic(fastq)
    // Post-trimming FastQC
    fastqc_post = fastqc_post(trimmed_fastq)
    // MultiQC report after trimming
    multiqc_input_post = fastqc_post.collect()
    multiqc_post = multiqc_post(multiqc_input_post)
    //--------------------------------------------
    // BWA_MEM mapping
    //--------------------------------------------
    aligned_bwa = bwa(trimmed_fastq, "${params.outdir}/sam_bwa")
    sorted_bwa = sort_bwa(aligned_bwa, "${params.outdir}/sorted_bwa")
    dedup_bwa = rm_dup_bwa(sorted_bwa, "${params.outdir}/dedup_bwa")
    map_bwa = mapped_bwa(dedup_bwa, "${params.outdir}/mapped_bwa")
    indexed_bwa = index_bwa(map_bwa, "${params.outdir}/mapped_bwa")
}
