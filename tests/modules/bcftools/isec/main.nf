#!/usr/bin/env nextflow



include { BCFTOOLS_ISEC } from '../../../../modules/bcftools/isec/main.nf'

workflow test_bcftools_isec {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true)],
              [ file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test2_vcf_gz_tbi'], checkIfExists: true)]
            ]

    BCFTOOLS_ISEC ( input )
}
