#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VSEARCH_SORT } from '../../../../../modules/nf-core/vsearch/sort/main.nf'

workflow test_vsearch_sort_size {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    VSEARCH_SORT ( input, "--sortbysize" )
}

workflow test_vsearch_sort_length {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    VSEARCH_SORT ( input, "--sortbylength" )
}
