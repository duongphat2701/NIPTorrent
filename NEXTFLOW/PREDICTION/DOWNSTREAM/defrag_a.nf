process defrag_a {
    input:
    path boy_dir
    path girl_dir
    path gcc
    path pickle
    
    output:
    path "*.tsv", emit: tsv

    // Publish results
    publishDir "${params.outdir}/defrag_a", mode: 'copy'

    script:
    """
    # Set up a safe temporary directory for R/Python tools
    export TMPDIR=\$(mktemp -d)
    echo "TMPDIR is set to \$TMPDIR"
    
    for file in ./*.gcc; do
        prefix=\$(basename \$file '.gcc')
        echo \$prefix >> prefix_sample.txt
    done

    ${projectDir}/bin/defrag_a.py prefix_sample.txt
    """
}