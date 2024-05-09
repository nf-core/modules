#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SHIFT } from '../../../../../modules/nf-core/bedtools/shift/main.nf'

workflow test_bedtools_shift {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]
    sizes = [[id:'sizes'],file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)]
    BEDTOOLS_SHIFT ( input, sizes )
}
