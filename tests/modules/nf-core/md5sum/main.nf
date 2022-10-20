#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MD5SUM } from '../../../../modules/nf-core/md5sum/main.nf'

workflow test_md5sum {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    MD5SUM ( input )
}
