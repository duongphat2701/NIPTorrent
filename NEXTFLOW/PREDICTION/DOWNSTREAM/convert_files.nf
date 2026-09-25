process convert_files {
    tag { sample_name }
    // Convert BAM files to .gcc, .pickle, .npz formats
    input:
    tuple val(sample_name), path(bam), path(bai)
    val outdir
    // Output channels
    output:
    tuple val(sample_name), path("*.pickle"), emit: pickle
    tuple val(sample_name), path("*.gcc"), emit: gcc
    tuple val(sample_name), path("*.npz"), emit: npz
    // Publish results
    publishDir "${params.outdir}/converted_files", mode: 'copy'
    // Script run
    script:
    """
    # Set up a safe temporary directory for R/Python tools
    export TMPDIR=\$(mktemp -d)
    ${projectDir}/bin/convert_files.py \\
        --binSizePickle ${params.binSize} \\
        --binSizeNpz ${params.binSize} \\
        ${bam}
    """
}
