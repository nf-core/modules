#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STITCH } from '../../../../modules/nf-core/stitch/main.nf'

workflow test_stitch {

    cramlist = Channel.fromPath(
        [
        params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram' ],
        params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'],
        ]
    )
    .map { it[-1] as String } // get only filename
    .collectFile( name: "cramlist.txt", newLine: true, sort: true )

    reads = Channel.of(
        [
            [ id:'test_reads' ], // meta map
            [
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'      ], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram_crai' ], checkIfExists: true),
            ],
            [
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'     ], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true),
            ],
        ]
    )
    .combine ( cramlist )

    reference = [
        [ id:'test_reference' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta']    , checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true),
    ]

    stitch_input = [
        [ id:'test_positions' ], // meta map
        file("../dbsnp_146.hg38.biallelic_snps.tsv", checkIfExists: true),
        [],
        [],
        "chr22",
        2,
        1
    ]

    STITCH ( stitch_input, reads, reference )
}
