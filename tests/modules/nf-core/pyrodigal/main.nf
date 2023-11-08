#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYRODIGAL } from '../../../../modules/nf-core/pyrodigal/main.nf'

workflow test_pyrodigal {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('/tmp/pytest_workflow_wb4hgjso/pyrodigal_test_pyrodigal/work/stage-0a3ce914-e002-407f-b594-c35213fecd3d/c9/9a152be072a23671b7a7961637229f/genome.fasta', checkIfExists: true)
    ]

    PYRODIGAL ( input )
}
