#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PEDDY } from '../../../modules/peddy/main.nf' addParams( options: [:] )

workflow test_peddy {

    input = [ [ id:'test', single_end:false ],
            file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)] // meta map
    vcf = file(params.test_data['homo_sapiens']['genome']['justhusky_minimal_vcf_gz'], checkIfExists: true)
    tbi = file(params.test_data['homo_sapiens']['genome']['justhusky_minimal_vcf_gz_tbi'], checkIfExists: true)

    PEDDY ( input , vcf , tbi )
}
