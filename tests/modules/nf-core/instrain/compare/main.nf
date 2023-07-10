#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INSTRAIN_COMPARE } from '../../../../../modules/phidepiper/instrain/compare/main.nf'
include { INSTRAIN_PROFILE } from '../../../../../modules/nf-core/instrain/profile/main.nf'
include { GUNZIP } from '../../../../../modules/nf-core/gunzip/main.nf'

workflow test_instrain_compare {
    
     input = Channel.of(
        [[ id:'test1', single_end:true ], file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)],
        [[ id:'test2', single_end:true ], file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)]
    )
    // input = Channel.of(
    //     [[ id:'test1', single_end:true ], file('./test_paired_end_sorted_bam1.sorted.bam')],
    //     [[ id:'test2', single_end:true ], file('./test_paired_end_sorted_bam2.sorted.bam')]
    // )


    genome_fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    

    INSTRAIN_PROFILE ( input , genome_fasta , [] , [] )
    INSTRAIN_PROFILE.out.profile.map{ [[id:"test"], it[1]] }.groupTuple().view()
    input.view()
    // map{ [[id:"test"], it[1] ] }.groupTuple().view()
    INSTRAIN_COMPARE ( INSTRAIN_PROFILE.out.profile.map{ [[id:"test"], it[1]] }.groupTuple() , input.map{ [[id:"test"], it[1]] }.groupTuple(), [] )

}
