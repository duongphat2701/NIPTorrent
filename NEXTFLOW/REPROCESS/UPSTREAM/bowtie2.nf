process bowtie2 {
    tag "${sample_name}"
    // Publish results
    publishDir "${params.outdir}/sam_bowtie2", mode: 'copy', overwrite: true

    input:
    tuple val(sample_name), path(read)

    output:
    tuple val(sample_name), path("${sample_name}.bam")

    script:
    """
    bowtie2 --very-sensitive --no-unal \
        -x ${params.ref_bowtie2} \
        -U ${read} \
        -S ${sample_name}.sam \
        -p ${task.cpus}
    samtools view -bS ${sample_name}.sam > ${sample_name}.bam
    """
}
