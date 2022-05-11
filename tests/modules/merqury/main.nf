#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MERYL_COUNT } from '../../../modules/meryl/count/main.nf'
include { MERQURY     } from '../../../modules/merqury/main.nf'

workflow test_merqury {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true)
    ]
    assembly = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    MERYL_COUNT ( input )
    MERQURY ( MERYL_COUNT.out.meryl_db.join( Channel.value( assembly ) ) )
}
