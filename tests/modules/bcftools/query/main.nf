#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_QUERY } from '../../../../modules/bcftools/query/main.nf' addParams( options: [:] )

workflow test_bcftools_query {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    BCFTOOLS_QUERY ( input )
}
