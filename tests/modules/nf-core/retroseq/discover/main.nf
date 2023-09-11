#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RETROSEQ_DISCOVER } from '../../../../../modules/nf-core/retroseq/discover/main.nf'

workflow test_retroseq_discover {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    RETROSEQ_DISCOVER ( input )
}
