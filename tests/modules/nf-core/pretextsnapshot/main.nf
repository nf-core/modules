#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRETEXTSNAPSHOT } from '../../../../modules/nf-core/pretextsnapshot/main.nf'

workflow test_pretextsnapshot {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PRETEXTSNAPSHOT_ALL ( input )
}
