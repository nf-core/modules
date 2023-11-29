#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_TBI } from '../../../../../modules/nf-core/tabix/bgziptabix/main.nf'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_CSI } from '../../../../../modules/nf-core/tabix/bgziptabix/main.nf'

workflow test_tabix_bgziptabix {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    TABIX_BGZIPTABIX_TBI ( input )
    TABIX_BGZIPTABIX_CSI ( input )

}
