#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQSCREEN_FASTQSCREEN } from '../../../../../modules/nf-core/fastqscreen/fastqscreen/main.nf'

workflow test_fastqscreen_fastqscreen {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    FASTQSCREEN_FASTQSCREEN ( input )
}
