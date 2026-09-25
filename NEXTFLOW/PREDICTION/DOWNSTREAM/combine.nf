process combine {
    publishDir "${outdir}", mode: 'copy'
    
    input:
    path tsvs
    val  output_name
    val  outdir

    output:
    path output_name, emit: summary

   
    // Script run
    script:
    """
    # Keep header from the first file
    head -n 1 ${tsvs[0]} > ${output_name}
    # Append all data rows (skip header in each file)
    tail -n +2 -q ${tsvs.join(' ')} >> ${output_name}
    """
}

