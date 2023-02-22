#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../../modules/nf-core/diamond/makedb/main.nf'
include { DIAMOND_BLASTP } from '../../../../../modules/nf-core/diamond/blastp/main.nf'

workflow test_diamond_blastp {

    db = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    out_ext = 'txt'
    blast_columns = 'qseqid qlen'

    DIAMOND_MAKEDB ( db )
    DIAMOND_BLASTP ( [ [id:'test'], fasta ], DIAMOND_MAKEDB.out.db, out_ext, blast_columns )
}

workflow test_diamond_blastp_daa {

    db = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    out_ext = 'daa'
    blast_columns = []

    DIAMOND_MAKEDB ( db )
    DIAMOND_BLASTP ( [ [id:'test'], fasta ], DIAMOND_MAKEDB.out.db, out_ext, blast_columns )
}
