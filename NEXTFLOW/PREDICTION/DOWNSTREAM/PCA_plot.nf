process pca_plot {
    // Generate PCA plots from converted files
    input:
    path gcc_files
    path pickle_files
    // Output channels
    output:
    path "pca.tsv", emit: tsv
    path "*pca.html"
    path "pca_plots"
    path "pca_bins"
    // Publish results
    publishDir "${params.outdir}/PCA", mode: 'copy'
    // Script run
    script:
    """
    # Set up a safe temporary directory for R/Python tools
    export TMPDIR=\$(mktemp -d)
    echo "TMPDIR is set to \$TMPDIR"
    
    for file in ./*.gcc; do
        prefix=\$(basename \$file '.gcc')
        echo \$prefix >> prefix_sample.txt
    done

    pca_plot.py prefix_sample.txt
    """
}