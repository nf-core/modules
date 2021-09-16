#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KLEBORATE } from '../../../modules/kleborate/main.nf' addParams( options: [:] )

workflow test_kleborate {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    KLEBORATE ( input )
}
