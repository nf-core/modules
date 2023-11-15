#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../modules/nf-core/diamond/makedb/main.nf'
include { EGGNOGMAPPER } from '../../../../modules/nf-core/eggnogmapper/main.nf'

workflow test_eggnogmapper {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    eggnog_db = file("https://github.com/nf-core/test-datasets/raw/eddf5b0e3336e0f93c81d4b4843b07257f9efaec/data/delete_me/eggnogmapper/eggnog.db", checkIfExists: true)
    eggnog_db.copyTo("${workDir}/tmp/eggnog.db")
    eggnog_data_dir = "${workDir}/tmp/"

    DIAMOND_MAKEDB ( [ [id:'test2'], fasta ] )
    EGGNOGMAPPER ( [ [id:'test'], fasta ], eggnog_db, eggnog_data_dir, DIAMOND_MAKEDB.out.db )
}
