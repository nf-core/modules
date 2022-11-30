#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR              } from '../../../../../modules/nf-core/untar/main.nf'
include { MMSEQS_CREATEINDEX } from '../../../../../modules/nf-core/mmseqs/createindex/main.nf'

workflow test_mmseqs_createindex {

    db = file('/home/kadomu/DELETE_ME/mmseqs_small_dbs/mmseqs.tar.gz', checkIfExists: true)

    ch_db = UNTAR ( [[:], db] ).untar
        .map { it[1] }
    MMSEQS_CREATEINDEX ( ch_db )
}
