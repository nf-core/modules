#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMBAMBA_MARKDUP } from '../../../../../modules/nf-core/sambamba/markdup/main.nf'

workflow test_sambamba_markdup {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SAMBAMBA_MARKDUP ( input )
}
