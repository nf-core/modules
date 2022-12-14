#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EMBOSS_GETORF } from '../../../../../modules/nf-core/emboss/getorf/main.nf'

workflow test_emboss_getorf {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)
    ]

    out_ext = 'orf'

    EMBOSS_GETORF ( input, out_ext )
}
