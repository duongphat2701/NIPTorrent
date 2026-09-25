process multiqc_pre {
    // Publish MultiQC report to output directory
    publishDir "${params.outdir}/multiqc_pre", mode: 'copy', overwrite: true
    // No input: MultiQC processes do not require input files directly
    input:
    path fastqc_pre, stageAs: "fastqc_pre/*"
    // Output: MultiQC report files
    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data
    // Script to run MultiQC
    script:
    """
    multiqc ${params.outdir}/fastqc_pre -o ./
    """
}
