#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FQ_FILTER } from '../../../../../modules/nf-core/fq/filter/main.nf'

workflow test_fq_filter {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    reads = Channel.of([]) // TODO: Populate

    FQ_FILTER ( input, empty )
}
