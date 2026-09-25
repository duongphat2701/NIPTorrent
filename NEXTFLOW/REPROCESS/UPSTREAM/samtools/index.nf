process index {
// Publish results
    publishDir "${outdir}", mode: 'copy', overwrite: true
// Intput: bam file
    input:
    tuple val(sample_name), path(bam)
    val outdir
// Output: bam file
    output:
    tuple val(sample_name), path("${sample_name}.bam.bai")
// Script: Run
    script:
    """
    samtools index -@ ${task.cpus} "${bam}"
    """
}
