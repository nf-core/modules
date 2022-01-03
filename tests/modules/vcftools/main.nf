#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFTOOLS as VCFTOOLS_BASE     } from '../../../modules/vcftools/main.nf'
include { VCFTOOLS as VCFTOOLS_OPTIONAL } from '../../../modules/vcftools/main.nf'

workflow test_vcftools_vcf_base {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
            ]

    VCFTOOLS_BASE ( input, [], [] )
}

workflow test_vcftools_vcfgz_base {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)
            ]

    VCFTOOLS_BASE ( input, [], [] )
}

workflow test_vcftools_vcf_optional {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
            ]
    bed   = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    VCFTOOLS_OPTIONAL ( input, bed, [] )
}

workflow test_vcftools_vcfgz_optional {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)
            ]
    bed   = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    VCFTOOLS_OPTIONAL ( input, bed, [] )
}
