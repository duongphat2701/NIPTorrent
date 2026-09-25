#!/usr/bin/env nextflow
include { convert_files } from "./DOWNSTREAM/convert_files.nf"
include { convert_gender } from "./DOWNSTREAM/convert_gender.nf"
include { create_ref } from "./DOWNSTREAM/create_ref.nf"
workflow {
    // -----------------------------
    // Load BAM and BAI
    // -----------------------------
    Channel
        .fromPath(params.input_csv)
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
        .flatMap { sample, bam, bai -> [bam, bai] }
        .collect()

    // ------------------------//
    //      Convert files      //
    // ------------------------//

    convert_files(bam_input)

    // ------------------------//
    //  Convert gender files   //
    // ------------------------//

    convert_gender(bam_input)

    // ------------------------//
    //   Create reference      //
    // ------------------------//

    create_ref(convert_files.out)
}
