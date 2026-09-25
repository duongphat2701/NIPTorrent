process create_ref {
    input:
    path gcc_files
    path pickle_files
    path npz_files
    path csv_table

    output:
    path "boydir", emit: boy_ref
    path "girldir", emit: girl_ref
    path "npz", emit: npz_ref
    path "wisecondorx_reference.npz", emit: npz
    // Publish results
    publishDir "${params.outdir}/reference", mode: 'copy'

    script:
    """
    for file in ./*.gcc; do
        prefix=\$(basename \$file '.gcc')
        echo \$prefix >> prefix_sample.txt
    done

    ${projectDir}/bin/create_ref.py prefix_sample.txt \
        --cpus ${task.cpus} \
        --ref-size ${params.refSize} \
        --wise-refsize ${params.wiseRefSize} \
        --binSizeNpz ${params.binSize}
    """
}
