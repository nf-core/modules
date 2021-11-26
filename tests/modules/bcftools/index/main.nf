#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_CSI } from '../../../../modules/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_TBI } from '../../../../modules/bcftools/index/main.nf'


workflow test_bcftools_index_csi {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ]

    BCFTOOLS_INDEX_CSI ( input )
}

workflow test_bcftools_index_tbi {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ]

    BCFTOOLS_INDEX_TBI ( input )
}
