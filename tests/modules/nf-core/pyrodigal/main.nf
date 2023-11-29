#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYRODIGAL } from '../../../../modules/nf-core/pyrodigal/main.nf'

workflow test_pyrodigal {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    PYRODIGAL ( input )
}
