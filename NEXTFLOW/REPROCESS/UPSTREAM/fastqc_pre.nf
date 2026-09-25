process fastqc_pre {
    tag "$sample_name"
    // Publish results
    publishDir "${params.outdir}/fastqc_pre", mode: 'copy', overwrite: true
    // Input: fastq file
    input:
    tuple val(sample_name), path(read)
    // Output: quality control report
    output:
    path('*_fastqc.{zip,html}')   
    // Script: Run
    script:
    """
    fastqc --threads ${task.cpus} ${read}
    """
}