#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTDENOVO } from '../../../../modules/nf-core/nextdenovo/main.nf'

workflow test_nextdenovo {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['nanopore']['test_fastq_gz'], checkIfExists: true)
    ]

    NEXTDENOVO ( input, "ont", "500k" )
}
