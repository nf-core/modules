#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAMSTATS_GENERALSTATS } from '../../../../../modules/nf-core/bamstats/generalstats/main.nf'

workflow test_bamstats_generalstats {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BAMSTATS_GENERALSTATS ( input )
}
