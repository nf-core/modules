#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_SORTSAM } from '../../../../modules/picard/sortsam/main.nf'

workflow test_picard_sortsam {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) ]
    sort_order = "queryname"

    PICARD_SORTSAM ( input, sort_order )
}
