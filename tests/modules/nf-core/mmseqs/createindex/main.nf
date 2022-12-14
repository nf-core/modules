#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR              } from '../../../../../modules/nf-core/untar/main.nf'
include { MMSEQS_CREATEINDEX } from '../../../../../modules/nf-core/mmseqs/createindex/main.nf'

workflow test_mmseqs_createindex {

    db = file(params.test_data['sarscov2']['genome']['mmseqs_tar_gz'], checkIfExists: true)

    ch_db = UNTAR ( [[:], db] ).untar
        .map { it[1] }
    MMSEQS_CREATEINDEX ( ch_db )
}
