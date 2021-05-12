#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

test_options = ['args': '--SCATTER_COUNT 6 --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW --UNIQUE true --SORT true']
include { GATK4_BEDTOINTERVALLIST } from '../../../../software/gatk4/bedtointervallist/main.nf' addParams( options: [:] )
include { GATK4_INTERVALLISTTOOLS as INTERVALLISTTOOLS } from '../../../../software/gatk4/intervallisttools/main.nf' addParams( options: test_options )

workflow test_gatk4_intervallisttools {

    input = [ [ id:'test' ], [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]]
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_BEDTOINTERVALLIST ( input, dict )
    INTERVALLISTTOOLS ( GATK4_BEDTOINTERVALLIST.out.interval_list )
}
