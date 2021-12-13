#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFANNO_LUA } from '../../../../modules/vcfanno/lua/main.nf'

workflow test_vcfanno_lua {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    VCFANNO_LUA ( input )
}
