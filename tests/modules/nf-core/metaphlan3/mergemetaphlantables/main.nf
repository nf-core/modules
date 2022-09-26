#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR                           } from '../../../../../modules/nf-core/untar/main.nf'
include { METAPHLAN3                      } from '../../../../../modules/nf-core/metaphlan3/metaphlan3/main.nf'
include { METAPHLAN3_MERGEMETAPHLANTABLES } from '../../../../../modules/nf-core/metaphlan3/mergemetaphlantables/main.nf'

workflow test_metaphlan3_mergemetaphlantables {

    input = Channel.of(
        [[ id:'test', single_end:true ], [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]],
        [[ id:'test2', single_end:true ], [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]]
    )

    db    = [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan_database.tar.gz', checkIfExists: true) ]

    UNTAR ( db )
    METAPHLAN3 ( input, UNTAR.out.untar.map{ it[1] } )
    METAPHLAN3_MERGEMETAPHLANTABLES ( METAPHLAN3.out.profile.map{ [[id:"test"], it[1]] }.groupTuple() )

}
