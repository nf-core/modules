#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPDISTS } from '../../../modules/snpdists/main.nf' addParams( options: [:] )

workflow test_snpdists {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    SNPDISTS ( input )
}
