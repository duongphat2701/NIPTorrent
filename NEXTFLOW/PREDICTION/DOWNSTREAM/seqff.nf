process seqff {
    tag "${sample}"
    // Calculate read counts from BAM files
    input:
    tuple val(sample), path(bam), path(bai)
    // Output channels
    output:
    path "${sample}.seqff.tsv", emit: tsv
    // Publish results
    publishDir "${params.outdir}/seqff", mode: 'copy'

    script:
    """
    # Set up a safe temporary directory for R/Python tools
    export TMPDIR=\$(mktemp -d)
    ${projectDir}/bin/seqff.py ${bam} ${sample}.seqff.tsv
    """
}

