#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SURVIVOR_MERGE } from '../../../../../modules/nf-core/survivor/merge/main.nf'

workflow test_survivor_merge {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
        ]
    ]

    SURVIVOR_MERGE ( input )
}
