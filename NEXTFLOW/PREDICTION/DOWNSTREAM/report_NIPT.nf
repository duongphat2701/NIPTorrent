process report {
    input:
    path read_count
    path pca
    path gender
    path seqff
    path abnormal
    path abnormal_path
    path plot_abnormal
    
    output:
    path "*.html", emit: html
    path "*.tsv", emit: tsv

    // Publish results
    publishDir "${params.outdir}/report", mode: 'copy'

    script:
    """
    report.py
    """
}