#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_FILTERBYNAME as BBMAP_FILTERBYNAME } from '../../../../../modules/nf-core/bbmap/filterbyname/main.nf'

workflow test_bbmap_filterbyname_paired_end {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    names = "ERR5069949.2151832,ERR5069949.576388,ERR5069949.501486"

    BBMAP_FILTERBYNAME ( input, names )
}

workflow test_bbmap_filterbyname_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    names = "ERR5069949.2151832,ERR5069949.576388,ERR5069949.501486"

    BBMAP_FILTERBYNAME ( input, names )
}
