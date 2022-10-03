#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYCOQC } from '../../../modules/pycoqc/main.nf'

workflow test_pycoqc {

    input = [ file(params.test_data['sarscov2']['nanopore']['test_sequencing_summary'], checkIfExists: true) ]

    PYCOQC ( input )
}
