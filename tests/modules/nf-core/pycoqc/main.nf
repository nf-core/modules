#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYCOQC } from '../../../../modules/nf-core/pycoqc/main.nf'

workflow test_pycoqc {

    def input = []
    input = [ [ id:'test' ], // meta map
               file(params.test_data['sarscov2']['nanopore']['test_sequencing_summary'], checkIfExists: true) ]

    PYCOQC ( input )
}
