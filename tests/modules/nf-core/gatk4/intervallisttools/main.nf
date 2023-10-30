#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_BEDTOINTERVALLIST } from '../../../../../modules/nf-core/gatk4/bedtointervallist/main.nf'
include { GATK4_INTERVALLISTTOOLS } from '../../../../../modules/nf-core/gatk4/intervallisttools/main.nf'

workflow test_gatk4_intervallisttools {

    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]
    dict = Channel.fromPath(file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true))
        .map{ dict -> [ [id:'test'], dict ] }

    GATK4_BEDTOINTERVALLIST ( input, dict )
    GATK4_INTERVALLISTTOOLS ( GATK4_BEDTOINTERVALLIST.out.interval_list )
}
