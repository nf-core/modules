#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../../modules/nf-core/diamond/makedb/main.nf'
include { DIAMOND_BLASTX } from '../../../../../modules/nf-core/diamond/blastx/main.nf'
include { MEGAN_DAA2INFO } from '../../../../../modules/nf-core/megan/daa2info/main.nf'

workflow test_megan_daa2info {

    db = [ file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['proteome_fasta'], checkIfExists: true) ]
    fasta = [ file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['genome_fasta'], checkIfExists: true) ]
    out_ext = 'daa'
    blast_columns = []
    megan_summary = true

    DIAMOND_MAKEDB ( db )
    DIAMOND_BLASTX ( [ [id:'test'], fasta ], DIAMOND_MAKEDB.out.db, out_ext, blast_columns )
    MEGAN_DAA2INFO ( DIAMOND_BLASTX.out.daa, megan_summary )
}
