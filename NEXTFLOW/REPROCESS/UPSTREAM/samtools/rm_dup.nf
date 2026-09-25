process rm_dup {
// Publish results
    publishDir "${outdir}", mode: 'copy', overwrite: true
// Remove duplicates from BAM files
    input:
    tuple val(sample_name), path(bam)
    val outdir
// Output: bam file
    output:
    tuple val(sample_name), path("${sample_name}_rmdup.bam")
// Script: Run
    script:
    """
    samtools markdup -@ ${task.cpus} -r "${bam}" "${sample_name}_rmdup.bam"
    rm "${bam}"
    """
}

