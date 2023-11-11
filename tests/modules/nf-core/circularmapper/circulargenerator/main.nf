#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CIRCULARMAPPER_CIRCULARGENERATOR } from '../../../../../modules/nf-core/circularmapper/circulargenerator/main.nf'

workflow test_circularmapper_circulargenerator {
    
    input = [
        [ id:'test', circularextension:500, circulartarget:'MT192765.1' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    CIRCULARMAPPER_CIRCULARGENERATOR ( input )
}
