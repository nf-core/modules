#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { WINDOWMASKER } from '../../../../modules/nf-core/windowmasker/main.nf'

workflow test_windowmasker {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    WINDOWMASKER ( [ [id:'test'], input ] )

}