process nipter {

    tag "${sample}"

    publishDir "${params.outdir}/NIPTeR", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith('_results.rds'))      "rds/$filename"
        else if (filename.endsWith('_summary_wide.tsv')) "abnormal_tables/$filename"
        else null
    }

    input:
    tuple val(sample), path(bam), path(bai), path(ref_nipter)

    output:
    path "${sample}.nipter_allchr_results.rds", emit: nipter_rds
    path "${sample}.nipter_allchr_summary_wide.tsv", emit: nipter_summary

    script:
    """
    RUN_NCV=TRUE Rscript ${projectDir}/bin/nipter1.R \
        ${ref_nipter} \
        ${bam} \
        ${sample}
    """
}
