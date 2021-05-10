#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

test_options = ['args': '--SCATTER_COUNT 6 --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW --UNIQUE true --SORT true']
include { GATK4_INTERVALLISTTOOLS as INTERVALLISTTOOLS_ALL } from '../../../../software/gatk4/intervallisttools/main.nf' addParams( options: test_options )

workflow test_gatk4_intervallisttools_all {

    input = [ [ id:'test.allparams' ], file('/home/praveen/Desktop/tools/downloads/human_gtf/human.exons.short.interval_list', checkIfExists: true) ]

    INTERVALLISTTOOLS_ALL ( input )
}
