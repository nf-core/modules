#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAXBIN2 } from '../../../modules/maxbin2/main.nf'

workflow test_maxbin2 {

    input = [ 
        [ id:'test1', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true),
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true),
        []
    ]

    MAXBIN2 ( input )
}
