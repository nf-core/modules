#!/usr/bin/env nextflow



include { PYCOQC } from '../../../modules/pycoqc/main.nf'

workflow test_pycoqc {

    input = [ file(params.test_data['sarscov2']['nanopore']['test_sequencing_summary'], checkIfExists: true) ]

    PYCOQC ( input )
}
