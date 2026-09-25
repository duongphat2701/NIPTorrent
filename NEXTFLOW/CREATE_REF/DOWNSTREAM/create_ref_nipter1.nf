process create_ref_nipter {

    tag "control_group"

    publishDir "${params.outdir}/NIPTeR/reference", mode: 'copy'

    input:
    path control_bam_list

    output:
    // These match nipter_pipeline.R --action build outputs
    path "${params.nipter_prefix ?: 'NIPTeR_ref'}.controlgroup.GC.rds",          emit: controlgroup_rds
    path "${params.nipter_prefix ?: 'NIPTeR_ref'}.controlgroup.diagnostics.rds", emit: diagnostics_rds

    // Optional outputs (only created if you enable flags in params)
    path "${params.nipter_prefix ?: 'NIPTeR_ref'}.ncv_templates.rds",            optional: true, emit: ncv_templates_rds
    path "${params.nipter_prefix ?: 'NIPTeR_ref'}.matchqc_stats.rds",            optional: true, emit: matchqc_stats_rds

    // Convenience: collect any RDS created by this step
    path "${params.nipter_prefix ?: 'NIPTeR_ref'}*.rds",                         emit: all_rds

    script:
    def prefix      = params.nipter_prefix ?: 'NIPTeR_ref'
    def gcMethod    = params.nipter_gc_method ?: 'LOESS'
    def ncvChroms   = params.nipter_ncv_chromosomes ?: '13,18,21'
    def mqSubset    = params.nipter_matchqc_subset_controls ?: 200

    def makeNcvFlag = (params.nipter_make_ncv_templates ?: false) \
        ? "--make_ncv_templates --ncv_chromosomes '${ncvChroms}'"
        : ""

    def mqStatsFlag = (params.nipter_build_matchqc_stats ?: false) \
        ? "--build_matchqc_stats --matchqc_subset_controls ${mqSubset}"
        : ""

    """
    Rscript ${projectDir}/bin/nipter_pipeline.R \
        --action build \
        --bam_list ${control_bam_list} \
        --outdir . \
        --out_prefix ${prefix} \
        --cpus ${task.cpus} \
        --separate_strands \
        --gc_method ${gcMethod} \
        --gc_include_xy \
        ${makeNcvFlag} \
        ${mqStatsFlag}
    """
}
