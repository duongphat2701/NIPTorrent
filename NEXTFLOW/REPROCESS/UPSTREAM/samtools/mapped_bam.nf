process mapped_bam {
// Publish results
    publishDir "${outdir}", mode: 'copy', overwrite: true
// Map BAM files to remove unmapped reads
    input:
    tuple val(sample_name), path(bam)
    val outdir
//  Output channels
    output:
    tuple val(sample_name), path("${sample_name}.bam")
// Script: Run
    script:
    """
    samtools view -@ ${task.cpus} -F 4 -b ${bam} > ${sample_name}.bam
    rm ${bam}
    """
}

