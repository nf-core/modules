#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_COMPLEMENT } from '../../../../modules/bedtools/complement/main.nf'

workflow test_bedtools_complement {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    BEDTOOLS_COMPLEMENT ( input, sizes )
}
