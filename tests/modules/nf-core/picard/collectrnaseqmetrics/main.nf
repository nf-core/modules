#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTRNASEQMETRICS } from '../../../../../modules/nf-core/picard/collectrnaseqmetrics/main.nf'

workflow test_picard_collectrnaseqmetrics {
    
    input = [
        [ id:'test', single_end:false,strandedness:'forward' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
    ]

    PICARD_COLLECTRNASEQMETRICS ( input, false, false, false )
}
