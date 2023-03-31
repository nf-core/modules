#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP    } from '../../../../modules/nf-core/gunzip/main.nf'
include { PYRODIGAL } from '../../../../modules/nf-core/pyrodigal/main.nf'

workflow test_pyrodigal {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    GUNZIP( input )
    PYRODIGAL ( GUNZIP.out.gunzip )
}
