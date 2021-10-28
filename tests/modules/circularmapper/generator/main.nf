#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CIRCULARMAPPER_GENERATOR } from '../../../../modules/circularmapper/generator/main.nf' addParams( options: [publish_dir:'circularmapper'] )

workflow test_circularmapper_generator {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    CIRCULARMAPPER_GENERATOR ( input )
}
