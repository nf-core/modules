#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CLAIR3 } from '../../../../modules/nf-core/clair3/main.nf'

workflow test_clair3 {

    bam_tuple_ch = Channel.of([ [ id:'test' ], // meta map
                               file(params.test_data['sarscov2']['nanopore']['test_sorted_bam'], checkIfExists: true),
                               file(params.test_data['sarscov2']['nanopore']['test_sorted_bam_bai'], checkIfExists: true),
                                ])

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    platform = "ont"
    clair3_model = "r941_prom_sup_g5014"

    CLAIR3 ( bam_tuple_ch, fasta, fai)
}
