process tmap {
    tag "$sample_name"
    // Publish results
    publishDir "${params.outdir}/sam_tmap", mode: 'copy', overwrite: true
    // Input: fastq file
    input:
    tuple val(sample_name), path(read)
    // Output: bam file
    output:
    tuple val(sample_name), path("${sample_name}.bam")
    // Script: Run
    script:
    """
    set -euo pipefail

    if [[ "${read}" == *.gz ]]; then
      zcat ${read} | \
        tmap mapall -f ${params.ref_tmap} -r - -i fastq -n ${task.cpus} -o 0 -v stage1 map4 | \
        samtools view -@ ${task.cpus} -bS - > ${sample_name}.bam
    else
      tmap mapall -f ${params.ref_tmap} -r ${read} -i fastq -n ${task.cpus} -o 0 -v stage1 map4 | \
        samtools view -@ ${task.cpus} -bS - > ${sample_name}.bam
    fi
    """
}