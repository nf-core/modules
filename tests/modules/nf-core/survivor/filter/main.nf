#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SURVIVOR_FILTER } from '../../../../../modules/nf-core/survivor/filter/main.nf'

workflow test_survivor_filter {
    
    input = [
        [ id:'test'], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
    ]

    SURVIVOR_FILTER ( 
        input,
        [],
        51,
        10001,
        0.01,
        10
    )
}
