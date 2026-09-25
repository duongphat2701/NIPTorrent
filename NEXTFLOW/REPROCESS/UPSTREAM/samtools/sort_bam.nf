process sort_bam {
    tag { sample_name }
// Publish results
    publishDir "${outdir}", mode: 'copy', overwrite: true
// Sort BAM files
    input:
    tuple val(sample_name), path(bam)
    val outdir
// Output: bam file
    output:
    tuple val(sample_name), path("${sample_name}_sorted.bam")
// Script: Run
    script:
    """
    samtools sort -@ ${task.cpus} -O bam "${bam}" -o "${sample_name}_sorted.bam"
    """
}
