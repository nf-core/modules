#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BACPHLIP } from '../../../../modules/nf-core/bacphlip/main.nf'

workflow test_bacphlip {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['scaffolds_fasta'], checkIfExists: true)
    ]

    BACPHLIP ( input )
}
