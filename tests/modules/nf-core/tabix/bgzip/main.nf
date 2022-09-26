#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_BGZIP } from '../../../../../modules/nf-core/tabix/bgzip/main.nf'

workflow test_tabix_bgzip_compress {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    TABIX_BGZIP ( input )
}

workflow test_tabix_bgzip_decompress {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed_gz'], checkIfExists: true) ]
            ]

    TABIX_BGZIP ( input )
}
