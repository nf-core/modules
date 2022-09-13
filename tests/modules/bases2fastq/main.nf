#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BASES2FASTQ } from '../../../modules/bases2fastq/main.nf'

workflow test_bases2fastq {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BASES2FASTQ ( input )
}
