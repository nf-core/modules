#!/usr/bin/env nextflow



include { BCFTOOLS_STATS } from '../../../../modules/bcftools/stats/main.nf'

workflow test_bcftools_stats {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    BCFTOOLS_STATS ( input )
}
