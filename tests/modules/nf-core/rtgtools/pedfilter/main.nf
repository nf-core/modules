#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RTGTOOLS_PEDFILTER } from '../../../../../modules/nf-core/rtgtools/pedfilter/main.nf'

workflow test_rtgtools_pedfilter_ped {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)
    ]

    RTGTOOLS_PEDFILTER ( input )
}

workflow test_rtgtools_pedfilter_vcf {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true)
    ]

    RTGTOOLS_PEDFILTER ( input )
}

workflow test_rtgtools_pedfilter_vcf_output {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)
    ]

    RTGTOOLS_PEDFILTER ( input )
}
