#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR                           } from '../../../../../modules/nf-core/untar/main.nf'
include { METAPHLAN_METAPHLAN             } from '../../../../../modules/nf-core/metaphlan/metaphlan/main.nf'
include { METAPHLAN_MERGEMETAPHLANTABLES  } from '../../../../../modules/nf-core/metaphlan/mergemetaphlantables/main.nf'

workflow test_metaphlan_mergemetaphlantables {

    input = Channel.of(
        [[ id:'test', single_end:true ], [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]],
        [[ id:'test2', single_end:true ], [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]]
    )

    db    = [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan4_database.tar.gz', checkIfExists: true) ]

    UNTAR ( db )
    METAPHLAN_METAPHLAN ( input, UNTAR.out.untar.map{ it[1] } )
    METAPHLAN_MERGEMETAPHLANTABLES ( METAPHLAN_METAPHLAN.out.profile.map{ [[id:"test"], it[1]] }.groupTuple() )

}
