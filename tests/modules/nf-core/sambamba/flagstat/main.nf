#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMBAMBA_FLAGSTAT } from '../../../../../modules/nf-core/sambamba/flagstat/main.nf'

workflow test_sambamba_flagstat {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SAMBAMBA_FLAGSTAT ( input )
}
