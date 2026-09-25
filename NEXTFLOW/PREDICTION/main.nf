#!/usr/bin/env nextflow
// Include downstream modules
include { convert_files }      from "./DOWNSTREAM/convert_files.nf"
include { read_count }         from "./DOWNSTREAM/read_count.nf"
include { pca_plot }           from "./DOWNSTREAM/PCA_plot.nf"
include { convert_gender}      from "./DOWNSTREAM/convert_gender.nf"
include { defrag_a }           from "./DOWNSTREAM/defrag_a.nf"
include { defrag_b }           from "./DOWNSTREAM/defrag_b.nf"
include { seqff }              from "./DOWNSTREAM/seqff.nf"
include { wisecondorX }        from "./DOWNSTREAM/wisecondorX.nf"
include { nipter }             from "./DOWNSTREAM/nipter.nf"
include { BLACKLIST_SIZE_MASK } from "./DOWNSTREAM/blacklist_size_mask.nf"
include { PRECI_ABNORMAL_WITH_BLACKLIST } from "./DOWNSTREAM/predict_abnormal_with_blacklist.nf"
include { PRECI_ABNORMAL_WITHOUT_BLACKLIST } from "./DOWNSTREAM/predict_abnormal_without_blacklist.nf"
include { combine as combine_defrag_b }         from "./DOWNSTREAM/combine.nf"
include { combine as combine_seqff }            from "./DOWNSTREAM/combine.nf"
include { combine as combine_wisecondorX }      from "./DOWNSTREAM/combine.nf"
include { combine as combine_nipter }           from "./DOWNSTREAM/combine.nf"
workflow {
outdir = params.outdir ?: "./results"
// -----------------------------
// Load BAM and BAI
// -----------------------------
    Channel .fromPath(params.input_csv)
            .splitCsv(header: true)
            .map { row ->
                def sample = row.sample
                def bam = file(row.bam)
                def bai = row.bai ? file(row.bai) : file("${row.bam}.bai")
                tuple(sample, bam, bai)
            }
            .set { bam_input }
    // Input: Channel: Collect bam and bai files
    bam_list = bam_input
        .flatMap { sample, bam, bai -> [ bam, bai ] }
        .collect()
// -----------------------------
// Define reference files
// -----------------------------
    Channel
        .fromPath("${params.reference_dir}/boydir")
        .set { boy_dir } 
    Channel
        .fromPath("${params.reference_dir}/girldir")
        .set { girl_dir }
    Channel
        .fromPath("${params.reference_dir}/wisecondorx_reference.npz")
        .set { ref_npz }
    // Channel
    //     .fromPath("${params.reference_dir}/NIPTeR_control.rds")
    //     .set { ref_nipter }
// -----------------------------
// DOWNSTREAM ANALYSES
// -----------------------------
// //
    // Convert bam file to .gcc, .pickle, .npz
    convert_files (bam_input, outdir)  
    converted_results_cluster = convert_files.out.gcc
            .combine(convert_files.out.pickle, by:0)
            .combine(convert_files.out.npz, by:0)
// //
    // READ COUNT
    read_count_results = read_count (bam_list)
// //
    // PCA PLOT
    pca_plot_results = pca_plot (
        converted_results_cluster.map { it[1] }.collect(),
        converted_results_cluster.map { it[2] }.collect()
    )
// //
    // GENDER PREDICTION
    convert_gender (bam_list, outdir)
    gender_csv = convert_gender.out.gender_prediction
// //
    // DEFRAG A
    defrag_a (
        boy_dir,
        girl_dir,
        converted_results_cluster.map { it[1] }.collect(),
        converted_results_cluster.map { it[2] }.collect()
    )
// //
    // DEFRAG B
    defrag_b_inputs = boy_dir
        .combine(girl_dir)
        .combine(converted_results_cluster)
        .map { boy, girl, sample, gcc, pickle, npz ->
            tuple(boy, girl, gcc, pickle)
        }
    defrag_b(
        defrag_b_inputs,
        gender_csv
    )
    combine_defrag_b(
        defrag_b.out.tsv.collect(),
        "defrag_b_summary.csv",
        "${params.outdir}/defrag_b"
    )
// //
    // SEQFF
    seqff(bam_input)
    combine_seqff(
        seqff.out.tsv.collect(),
        "seqff_summary.csv",
        "${params.outdir}/seqff"
    )
// //
    // WISECONDORX
    wisecondorX(
        convert_files.out.npz,
        ref_npz
    )
    combine_wisecondorX(
        wisecondorX.out.abnormal_table.collect(),
        "wisecondorX_summary.csv",
        "${params.outdir}/wisecondorX"
    )
// //
    // WISECONDORX without MASK
    predict_abnormal_without_blacklist_results = PRECI_ABNORMAL_WITHOUT_BLACKLIST (
        convert_files.out.npz,
        ref_npz,
        "${outdir}/abnormal_without_blacklist_prediction"
    )
    collected_z_bins = predict_abnormal_without_blacklist_results.z_bin.collect()
    // WISECONDORX with MASK
    blacklist_size_mask = BLACKLIST_SIZE_MASK (
        collected_z_bins,
        params.basic_blacklist,
        outdir
    ) 
    merged_input_blacklist = convert_files.out.npz.join(blacklist_size_mask.blacklist_file)
    predict_abnormal_with_blacklist_results = PRECI_ABNORMAL_WITH_BLACKLIST (
        merged_input_blacklist,
        ref_npz,
        "${outdir}/abnormal_size_mask_blacklist"
    )
    // NIPTeR
    // samples_with_ref = bam_input.combine(ref_nipter)
    // nipter(samples_with_ref)
    // combine_nipter(
    //     nipter.out.nipter_summary.collect(),
    //     "nipter_summary.csv",
    //     "${params.outdir}/NIPTeR"
    // )
}