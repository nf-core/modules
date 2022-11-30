#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../../modules/nf-core/diamond/makedb/main.nf'
include { DIAMOND_BLASTX } from '../../../../../modules/nf-core/diamond/blastx/main.nf'

workflow test_diamond_blastx {

    db = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    fasta = [ file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]
    out_ext = 'tfdfdt'  // Nonsense file extension to check default case.
    blast_columns = 'qseqid qlen'

    DIAMOND_MAKEDB ( db )
    DIAMOND_BLASTX ( [ [id:'test'], fasta ], DIAMOND_MAKEDB.out.db, out_ext, blast_columns )
}

workflow test_diamond_blastx_daa {

    db = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    fasta = [ file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]
    out_ext = 'daa'
    blast_columns = []

    DIAMOND_MAKEDB ( db )
    DIAMOND_BLASTX ( [ [id:'test'], fasta ], DIAMOND_MAKEDB.out.db, out_ext, blast_columns )
}
