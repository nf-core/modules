#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTHSMETRICS } from '../../../../modules/picard/collecthsmetrics/main.nf' addParams( options: [:] )

workflow test_picard_collecthsmetrics {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    PICARD_COLLECTHSMETRICS ( input )
}
