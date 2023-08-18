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
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta']    , checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true),
    ]

    stitch_input = [
        [ id:'test_positions' ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/sequence/dbsnp_138.hg38.first_10_biallelic_sites.tsv", checkIfExists: true),
        [],
        [],
        "chr21",
        2,
        1
    ]

    STITCH ( stitch_input, reads, reference )
}
