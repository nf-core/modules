#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVPYTOR_PARTITION } from '../../../../modules/cnvpytor/partition/main.nf'

workflow test_cnvpytor_partition {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_pytor'], checkIfExists: true)
    ]

    CNVPYTOR_PARTITION ( input )
}
