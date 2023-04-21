#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQSCREEN_BUILD } from '../../../../../modules/nf-core/fastqscreen/build/main.nf'

workflow test_fastqscreen_build {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    FASTQSCREEN_BUILD ( input )
}
