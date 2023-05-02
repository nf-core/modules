#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTDENOVO } from '../../../../modules/nf-core/nextdenovo/main.nf'

workflow test_nextdenovo {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)
    ]

    NEXTDENOVO ( input, "hifi", "0.2k" )
}
