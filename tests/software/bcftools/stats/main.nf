#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_STATS } from '../../../../software/bcftools/stats/main.nf' addParams( options: [:] )

workflow test_bcftools_stats {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    BCFTOOLS_STATS ( input )
}
