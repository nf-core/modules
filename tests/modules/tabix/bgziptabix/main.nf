#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_BGZIPTABIX } from '../../../../modules/tabix/bgziptabix/main.nf' addParams( options: ['args2': '-p vcf'] )

workflow test_tabix_bgziptabix {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    TABIX_BGZIPTABIX ( input )
}
