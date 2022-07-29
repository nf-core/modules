#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_INTERVALLISTTOBED } from '../../../../modules/gatk4/intervallisttobed/main.nf'

workflow test_gatk4_intervallisttobed {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list'], checkIfExists: true)
    ]

    GATK4_INTERVALLISTTOBED ( input )
}
