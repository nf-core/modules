#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SLOP } from '../../../../software/bedtools/slop/main.nf' addParams( options: [args: '-l 15 -r 30', suffix: '_out'] )

workflow test_bedtools_slop {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)
    
    BEDTOOLS_SLOP ( input, sizes )
}

