#!/usr/bin/env nextflow



include { TABIX_BGZIPTABIX } from '../../../../modules/tabix/bgziptabix/main.nf'

workflow test_tabix_bgziptabix {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    TABIX_BGZIPTABIX ( input )
}
