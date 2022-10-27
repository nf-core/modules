#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { TABIX_BGZIP } from "$moduleDir/modules/nf-core/tabix/bgzip/main.nf"
include { TABIX_BGZIP as TABIX_BGZIP_WITH_GZI } from "$moduleDir/modules/nf-core/tabix/bgzip/main.nf"

workflow test_tabix_bgzip_compress {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    TABIX_BGZIP ( input )
}

workflow test_tabix_bgzip_compress_gzi {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    TABIX_BGZIP_WITH_GZI ( input )
}

workflow test_tabix_bgzip_decompress {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed_gz'], checkIfExists: true) ]
            ]

    TABIX_BGZIP ( input )
}
