// process create_ref_nipter {
//     tag "control_group"

//     publishDir "${params.outdir}/NIPTeR", mode: 'copy'

//     input:
//     path control_bam_list  // this is the 'control_bams.txt' file from above

//     output:
//     path "*.rds", emit: nipter_ref
//     script:
//     """
//     Rscript ${projectDir}/bin/ref_nipter.R \
//         --bam_list ${control_bam_list} \
//         --out_rds NIPTeR_control.rds
//     """
// }
process create_ref_nipter {

    tag "control_group"

    publishDir "${params.outdir}/NIPTeR/reference", mode: 'copy'

    input:
    path bam_files
    path bam_index

    output:
    path "NIPTeR_ref.controlgroup.GC.rds",           emit: controlgroup_rds
    path "NIPTeR_ref.controlgroup.diagnostics.rds",  emit: diagnostics_rds
    path "NIPTeR_ref.ncv_templates.rds",             optional: true, emit: ncv_templates_rds
    path "NIPTeR_ref.matchqc_stats.rds",             optional: true, emit: matchqc_stats_rds
    
    script:
    // Enable Z + NCV + RBZ + chrX/Y Z + QC
    def make_ncv = params.containsKey('nipter_make_ncv_templates') ? params.nipter_make_ncv_templates : false
    def make_mq  = params.containsKey('nipter_build_matchqc_stats') ? params.nipter_build_matchqc_stats : false
    def ncv_chrs  = params.nipter_ncv_chromosomes ?: "13,18,21"
    def mq_subset = params.nipter_matchqc_subset_controls ?: 200
    def ncv_flag = make_ncv ? "--make_ncv_templates --ncv_chromosomes '${ncv_chrs}'" : ""
    def mq_flag  = make_mq  ? "--build_matchqc_stats --matchqc_subset_controls ${mq_subset}" : ""

    """
    ls -1 *.bam | sort > control_bams.txt

    Rscript ${projectDir}/bin/ref_nipter1.R \
      --bam_list control_bams.txt \
      --outdir . \
      --out_prefix NIPTeR_ref \
      --cpus ${task.cpus} \
      --gc_method LOESS \
      --gc_include_xy \
      ${ncv_flag} \
      ${mq_flag}
    """
}
