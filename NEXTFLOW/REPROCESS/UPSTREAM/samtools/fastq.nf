process convert_fastq {
    tag "$sample_name"
// Publish results
    publishDir("${params.outdir}/fastq", mode: 'copy', overwrite: true)
// Input: bam file
    input:
    tuple val(sample_name), path(read)
// Output: fastq file
    output:
    tuple val(sample_name), path("${sample_name}.fastq.gz")
// Script: Run
    script:
    """
    samtools fastq -@ ${task.cpus} "${read}" | gzip > "${sample_name}.fastq.gz"
    """
}
