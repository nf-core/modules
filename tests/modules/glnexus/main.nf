#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLNEXUS } from '../../../modules/glnexus/main.nf'

workflow test_glnexus {
    input = [
        [ id:'test' ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf_gz'], checkIfExists: true)
        ]
    ]

    GLNEXUS ( input )
}
