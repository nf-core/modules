#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_CONCAT } from '../../../../modules/bcftools/concat/main.nf' addParams( options: ['args': '--force-samples --no-version'] )

workflow test_bcftools_concat {
    
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true) ]
            ]


    BCFTOOLS_CONCAT ( input )
}
