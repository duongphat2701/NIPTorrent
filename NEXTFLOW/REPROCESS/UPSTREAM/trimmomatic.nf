process trimmomatic {
    tag "$sample_name"
    // Publish results
    publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true
    // Input: fastq file
    input:
    tuple val(sample_name), path(read)
    // Output: trimmed fastq file
    output:
    tuple val(sample_name), path("${sample_name}.fastq.gz")
    // Script: Run
    script:
    """
    trimmomatic SE -threads ${task.cpus} \
    "${read}" "${sample_name}_trim.fastq" \
    ${params.trimmomatic_options}

    gzip -c "${sample_name}_trim.fastq" > "${sample_name}.fastq.gz"
    rm "${sample_name}_trim.fastq"
    """
}
   // trimmomatic SE -threads ${task.cpus} \
    // ${convert_fastq} ${sample_name}.fastq \
    // SLIDINGWINDOW:50:15 MINLEN:50
    // rm ${convert_fastq}