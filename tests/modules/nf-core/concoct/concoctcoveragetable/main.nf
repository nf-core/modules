#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CONCOCT_CUTUPFASTA } from '../../../../../modules/nf-core/concoct/cutupfasta/main.nf'
include { CONCOCT_CONCOCTCOVERAGETABLE } from '../../../../../modules/nf-core/concoct/concoctcoveragetable/main.nf'

workflow test_concoct_concoctcoveragetable {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    produce_bedfile = true

    CONCOCT_CUTUPFASTA ( input, produce_bedfile )

    input2 = Channel.fromList([
        [
            [ id:'test', single_end:false ], // meta map
            [
                file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ],
            [
                file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]
        ]
    ])

    ch_concoctcoveragetable_input = CONCOCT_CUTUPFASTA.out.bed.join( input2 )

    CONCOCT_CONCOCTCOVERAGETABLE ( ch_concoctcoveragetable_input )
}
