process wisecondorX {

    tag { sample_name }

    publishDir "${params.outdir}/wisecondorX", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith('_aberrations.bed'))      "aberrations/$filename"
        else if (filename.endsWith('_bins.bed'))        "z_bins/$filename"
        else if (filename.endsWith('_segments.bed'))    "z_segments/$filename"
        else if (filename.endsWith('_statistics.txt'))  "statistics/$filename"
        else if (filename.endsWith('_abn.tsv'))         "abnormal_tables/$filename"
        else if (filename.endsWith('.plots'))           "plots/$filename"
        else if (filename.endsWith('_plot_tmp.json'))   "jsons/$filename"
        else null
    }

    input:
    tuple val(sample_name), path(npz_file)
    each path(ref_npz)

    output:
    path("${sample_name}_aberrations.bed"),                      emit: aberrations
    tuple val(sample_name), path("${sample_name}_bins.bed"),     emit: z_bins
    tuple val(sample_name), path("${sample_name}_segments.bed"), emit: z_segments
    path("${sample_name}_statistics.txt"),                       emit: statistics
    path("${sample_name}_abn.tsv"),                              emit: abnormal_table
    path("${sample_name}.plots"),                                emit: plot
    path("${sample_name}_plot_tmp.json"),                        emit: json

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

    python3 ${projectDir}/bin/wisecondorX.py "\$NPZ" --ref "${ref_npz}"
    """
}