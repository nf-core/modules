#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { JASMINESV } from '../../../../modules/nf-core/jasminesv/main.nf'

workflow test_jasminesv {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
        ],
        []
    ]

    JASMINESV ( input )
}
