#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CIRCULARMAPPER_GENERATOR } from '../../../../modules/circularmapper/generator/main.nf' addParams( options: [:] )

workflow test_circularmapper_generator {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    CIRCULARMAPPER_GENERATOR ( input )
}
