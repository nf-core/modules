#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_GATHERPILEUPSUMMARIES } from '../../../../modules/gatk4/gatherpileupsummaries/main.nf'

workflow test_gatk4_gatherpileupsummaries {

    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test_pileups_table'], checkIfExists: true)]
        //file(params.test_data['homo_sapiens']['illumina']['test_pileups_table'], checkIfExists: true)]
    ]

    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    GATK4_GATHERPILEUPSUMMARIES ( input, dict )
}
