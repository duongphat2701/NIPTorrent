process convert_gender {
    // Predict gender from BAM files
    input:
    path bam_list, stageAs:"bam_list/*"
    val outdir
    // Output channels
    output:
    path "gender_prediction_Yfrac.tsv", emit: yfrac
    path "gender_prediction.csv", emit: gender_prediction
    path "threshold.txt", emit: txt, optional: true
    // Publish results
    publishDir "${params.outdir}/gender_prediction", mode: 'copy'
    // Script run
    script:
    """
    for file in bam_list/*.bam; do
        echo \$file >> bam_file.txt
    done

    ${projectDir}/bin/convert_gender.py bam_file.txt -t ${params.gender_threshold}
    echo "${params.gender_threshold}" > threshold.txt 
    """
}
