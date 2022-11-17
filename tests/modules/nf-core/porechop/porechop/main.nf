#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PORECHOP_PORECHOP } from '../../../../../modules/nf-core/porechop/porechop/main.nf'

workflow test_porechop_porechop {

    input = [ [ id:'test', single_end:true ], // meta map
                file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]

    PORECHOP_PORECHOP ( input )
}
