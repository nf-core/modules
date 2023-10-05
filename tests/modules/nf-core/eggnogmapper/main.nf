#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../modules/nf-core/diamond/makedb/main.nf'
include { EGGNOGMAPPER } from '../../../../modules/nf-core/eggnogmapper/main.nf'

workflow test_eggnogmapper {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    eggnog_db = file("tests/modules/nf-core/eggnogmapper/data/fixtures/eggnog.db", checkIfExists: true)
    eggnog_data_dir = eggnog_db.parent

    DIAMOND_MAKEDB ( fasta )
    EGGNOGMAPPER ( [ [id:'test'], fasta ], eggnog_db, eggnog_data_dir, DIAMOND_MAKEDB.out.db )
}
