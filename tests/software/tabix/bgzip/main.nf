#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_BGZIP } from '../../../../software/tabix/bgzip/main.nf' addParams( options: [:] )

workflow test_tabix_bgzip {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    TABIX_BGZIP ( input )
}
