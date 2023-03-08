#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREYJA_VARIANTS } from '../../../../../modules/nf-core/freyja/variants/main.nf'

workflow test_freyja_variants {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    FREYJA_VARIANTS ( input )
}
