#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYCOQC } from '../../../software/pycoqc/main.nf' addParams ( options: ['args' : '--min_pass_qual 0'] )

workflow test_pycoqc {

    input = [ file(params.test_data['sarscov2']['nanopore']['test_sequencing_summary'], checkIfExists: true) ]

    PYCOQC ( input )
}
