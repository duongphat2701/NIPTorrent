#!/usr/bin/env nextflow
include { convert_files } from "./DOWNSTREAM/convert_files.nf"
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
        .flatMap { sample, bam, bai -> [ bam, bai ] }
        .collect()
    // Input: Channel: Get bam directories   
    bamdir = bam_input
        .map { sample, bam, bai -> bam.parent.toString() }
        .unique()

    // -----------------------------
    // Create sanefalcon reference
    // -----------------------------
    
        create_ref_sanefalcon( bamdir.first() )

    /// -----------------------------
    // Create NIPTeR control group
    // -----------------------------

        control_bam_list = bam_input
            .map { sample, bam, bai -> bam }
            .collectFile(name: 'control_bams.txt') { it.toString() + '\n' }

        create_ref_nipter(control_bam_list)

    // -----------------------------
    // Create sanefalcon reference
    // -----------------------------

        create_ref_sanefalcon(bam_list)

    // ------------------------//
    // Create gender reference //
    // ------------------------//

    convert_files(bam_input)

}
