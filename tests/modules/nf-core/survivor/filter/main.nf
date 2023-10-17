#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SURVIVOR_FILTER  as SURVIVOR_FILTER_NO_BED} from '../../../../../modules/nf-core/survivor/filter/main.nf'
include { SURVIVOR_FILTER  as SURVIVOR_FILTER_BED} from '../../../../../modules/nf-core/survivor/filter/main.nf'

workflow test_survivor_filter_no_bed {
    
    input_no_bed = [
        [ id:'test'], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
            []
    ]

    SURVIVOR_FILTER_NO_BED ( 
        input_no_bed,
        51,
        10001,
        0.01,
        10
    )

}

workflow test_survivor_filter {

    input_bed = [
        [ id:'test'], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
        ]

    SURVIVOR_FILTER_BED ( 
        input_bed,
        51,
        10001,
        0.01,
        10
    )

}
