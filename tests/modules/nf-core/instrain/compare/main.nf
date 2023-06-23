#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP            } from '../../../../../modules/nf-core/gunzip/main.nf'
include { INSTRAIN_PROFILE  } from '../../../../../modules/nf-core/instrain/profile/main.nf'
include { INSTRAIN_COMPARE  } from '../../../../../modules/nf-core/instrain/compare/main.nf'

workflow test_instrain_compare {

    input = Channel.of(
        [[ id:'test1', single_end:true ], file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)],
        [[ id:'test2', single_end:true ], file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)],
    )

    combined_bams = [
        [ id: 'test', single_end:false ],
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
        [ file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true) ],
    ]

    genome_fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    INSTRAIN_PROFILE ( input , genome_fasta , [] , [] )
    INSTRAIN_PROFILE.out.profile.map{ [[id:"test"], it[1]] }.groupTuple().view()
    input.map{ [[id:"test"], it[1] ] }.groupTuple().view()
    INSTRAIN_COMPARE ( INSTRAIN_PROFILE.out.profile.map{ [[id:"test"], it[1]] }.groupTuple() , input.map{ [[id:"test"], it[1]] }.groupTuple() , [] )
}
