#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_INDEX } from '../../../../modules/bcftools/index/main.nf' addParams( options: [:] )

workflow test_bcftools_index {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ]

    BCFTOOLS_INDEX ( input )
}
