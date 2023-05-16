#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_BGZIPTABIX as BGZIPTABIX_TBI } from '../../../../../modules/nf-core/tabix/bgziptabix/main.nf'
include { TABIX_BGZIPTABIX as BGZIPTABIX_CSI } from '../../../../../modules/nf-core/tabix/bgziptabix/main.nf'

workflow test_tabix_bgziptabix {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    BGZIPTABIX_TBI ( input )
    BGZIPTABIX_CSI ( input )

}
