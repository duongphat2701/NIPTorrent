process read_count {
    // Calculate read counts from BAM files
    input:
    path bam_list, stageAs:"bam_list/*"
    // Output channels
    output:
    path "*.tsv", emit: tsv
    // Publish results
    publishDir "${params.outdir}/read_count", mode: 'copy'
    // Script run
    script:
    """
    # Set up a safe temporary directory for R/Python tools
    export TMPDIR=\$(mktemp -d)
    echo "TMPDIR is set to \$TMPDIR"
    
    for file in bam_list/*.bam; do
        echo \$file >> bam_file.txt
    done

    ${projectDir}/bin/read_count.py bam_file.txt
    """
}