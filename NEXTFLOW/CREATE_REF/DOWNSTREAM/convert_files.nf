process convert_files {
    // Publish results
    publishDir "${params.outdir}/converted_files", mode: 'copy'
    
    tag { sample }
    // Convert BAM files to .gcc, .pickle, .npz formats
    input:
    tuple val(sample), path(bam), path(bai)
    // Output channels
    output:
    tuple val(sample), path("*.pickle"), emit: pickle
    tuple val(sample), path("*.gcc"), emit: gcc
    tuple val(sample), path("*.npz"), emit: npz
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
