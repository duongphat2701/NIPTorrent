process defrag_b {
    input:
    tuple path(boy_dir), path(girl_dir), val(gcc), val(pickle)
    path (gender_csv)
    
    output:
    path "*.tsv", emit: tsv
    // Publish results
    publishDir "${params.outdir}/defrag_b", mode: 'copy'

    script:
    """
    ${projectDir}/bin/defrag_b.py ${pickle} ${girl_dir} ${pickle.baseName}.defrag_b.tsv ${gender_csv}
    """
}