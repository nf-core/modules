#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

test_options = ['args': '-w 50 ']
include { BEDTOOLS_MAKEWINDOWS } from '../../../../modules/bedtools/makewindows/main.nf' addParams( options: test_options )

workflow test_bedtools_makewindows {
    
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)]

    BEDTOOLS_MAKEWINDOWS ( input, true )
}
