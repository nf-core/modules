#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEDUP } from '../../../../modules/nf-core/dedup/main.nf'

workflow test_dedup {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    DEDUP ( input )
}
