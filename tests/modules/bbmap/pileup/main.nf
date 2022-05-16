#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_PILEUP } from '../../../../modules/bbmap/pileup/main.nf'

workflow test_bbmap_pileup {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BBMAP_PILEUP ( input )
}
