#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_QUERY } from '../../../../modules/bcftools/query/main.nf' addParams( options: [:] )

workflow test_bcftools_query {
    
    input = [[ id:'test2', single_end:false ], // meta map
             [ file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)]]


    BCFTOOLS_QUERY ( input )
}
