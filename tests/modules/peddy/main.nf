#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PEDDY } from '../../../modules/peddy/main.nf' addParams( options: [:] )

workflow test_peddy {

    input = [
        [ id:'test', single_end:false ],
        file(params.test_data['homo_sapiens']['genome']['justhusky_minimal_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['justhusky_minimal_vcf_gz_tbi'], checkIfExists: true)
    ]
    ped = file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)

    PEDDY ( input , ped )
}
