process PRECI_ABNORMAL_WITHOUT_BLACKLIST {
    tag "${sample_name}"
    publishDir "${outdir}/abnormal", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith('_aberrations.bed')) "aberrations/$filename"
        else if (filename.endsWith('_bins.bed')) "bins/$filename"
        else if (filename.endsWith('_segments.bed')) "segments/$filename"
        else if (filename.endsWith('_statistics.txt')) "statistics/$filename"
        else if (filename.endsWith('.plots')) "plots/$filename"
        else null
    }
    input:
    tuple val(sample_name), path (npz_file)
    path ref_npz
    val outdir

    output:
    path ("${sample_name}_aberrations.bed"), emit: abn_seg
    tuple val(sample_name), path ("${sample_name}_bins.bed"), emit: z_bin
    tuple val(sample_name), path ("${sample_name}_segments.bed"), emit: z_seg
    path ("${sample_name}_statistics.txt"), emit: z_chr
    path ("${sample_name}_abn.tsv"), emit: tsv
    path ("${sample_name}.plots"), emit: plot
    script:
    """
    TARGET="${sample_name}.npz"
    if [ "\$(basename "${npz_file}")" = "\$TARGET" ]; then
        NPZ="${npz_file}"
    else
        SRC="\$(readlink -f "${npz_file}")"
        ln -sfn "\$SRC" "\$TARGET"
        NPZ="\$TARGET"
    fi
    
    ${projectDir}/bin/preci_abnormal.py $npz_file --ref ${ref_npz}
    """
}
