#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TMB } from '../../../modules/tmb/main.nf'

workflow test_tmb {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true)
    ]

    TMB ( input, "", "./config/mutect2",[] )
}
