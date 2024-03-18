#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQWISH_INDUCE } from '../../../../../modules/nf-core/seqwish/induce/main.nf'

workflow test_seqwish_induce_transcriptome {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['transcriptome_paf'], checkIfExists: true)],
              [ file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]
            ]

    SEQWISH_INDUCE ( input )
}

workflow test_seqwish_induce_pangenome {
    input = [ [ id:'test' ], // meta map
            [ file(params.test_data['homo_sapiens']['pangenome']['pangenome_paf'], checkIfExists: true)],
            [ file(params.test_data['homo_sapiens']['pangenome']['pangenome_fa'], checkIfExists: true) ]
        ]

    SEQWISH_INDUCE ( input )
}
