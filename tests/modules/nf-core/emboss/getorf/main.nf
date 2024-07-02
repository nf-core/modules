#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EMBOSS_GETORF } from '../../../../../modules/nf-core/emboss/getorf/main.nf'

workflow test_emboss_getorf {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['transcriptome_fasta'], checkIfExists: true)
    ]

    out_ext = 'fasta'

    EMBOSS_GETORF ( input, out_ext )
}
