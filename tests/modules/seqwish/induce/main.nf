#!/usr/bin/env nextflow



include { SEQWISH_INDUCE } from '../../../../modules/seqwish/induce/main.nf'

workflow test_seqwish_induce {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['transcriptome_paf'], checkIfExists: true)],
              [ file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]
            ]

    SEQWISH_INDUCE ( input )
}
