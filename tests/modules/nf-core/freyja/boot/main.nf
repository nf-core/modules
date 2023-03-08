#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREYJA_BOOT } from '../../../../../modules/nf-core/freyja/boot/main.nf'

workflow test_freyja_boot {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    FREYJA_BOOT ( input )
}
