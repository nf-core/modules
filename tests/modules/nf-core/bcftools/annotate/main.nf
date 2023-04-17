#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_VCF } from '../../../../../modules/nf-core/bcftools/annotate/main.nf'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_BCF } from '../../../../../modules/nf-core/bcftools/annotate/main.nf'

workflow test_bcftools_annotate_out_vcf {

    input = [ [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true),
        [],[]
    ]

    BCFTOOLS_ANNOTATE_VCF ( input, [] )
}

workflow test_bcftools_annotate_out_bcf {

    input = [ [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_bcf'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true),
        [],[]
    ]

    BCFTOOLS_ANNOTATE_BCF ( input, [] )
}
