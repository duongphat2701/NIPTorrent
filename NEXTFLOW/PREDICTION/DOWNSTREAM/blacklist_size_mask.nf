process BLACKLIST_SIZE_MASK {
    tag "generate_blacklist"
    publishDir "${outdir}/blacklist_size_mask", mode: 'copy'

    input:
    tuple val(sample_name), path(input_bins)
    path basic_blacklist
    val outdir

    output:
    tuple val(sample_name), path("*_size_mask_blacklist.bed"), emit: blacklist_file

    script:
    """
    mkdir -p input_${sample_name}
    cp ${input_bins} input_${sample_name}/

    ${projectDir}/bin/generate_size_mask_blacklist.py \\
        --input-folder "input_${sample_name}" \\
        --output-folder . \\
        --basic-blacklist "$basic_blacklist"
    """
}
