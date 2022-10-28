#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CONCOCT_CUTUPFASTA           } from '../../../../../modules/nf-core/concoct/cutupfasta/main.nf'
include { CONCOCT_CONCOCTCOVERAGETABLE } from '../../../../../modules/nf-core/concoct/concoctcoveragetable/main.nf'
include { CONCOCT_CONCOCT              } from '../../../../../modules/nf-core/concoct/concoct/main.nf'
include { CONCOCT_MERGECUTUPCLUSTERING } from '../../../../../modules/nf-core/concoct/mergecutupclustering/main.nf'
include { CONCOCT_EXTRACTFASTABINS     } from '../../../../../modules/nf-core/concoct/extractfastabins/main.nf'

workflow test_concoct_extractfastabins {

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

    ch_concoctconcoct_input = CONCOCT_CONCOCTCOVERAGETABLE.out.tsv
                                .join(CONCOCT_CUTUPFASTA.out.fasta)

    CONCOCT_CONCOCT( ch_concoctconcoct_input )

    CONCOCT_MERGECUTUPCLUSTERING ( CONCOCT_CONCOCT.out.clustering_csv )

    ch_input_for_extractfastabins = Channel.of(input).join(CONCOCT_MERGECUTUPCLUSTERING.out.csv)

    CONCOCT_EXTRACTFASTABINS ( ch_input_for_extractfastabins )

    CONCOCT_EXTRACTFASTABINS.out.fasta
}
