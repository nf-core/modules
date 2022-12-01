#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR                  } from '../../../../../modules/nf-core/untar/main.nf'
include { MMSEQS_TSV2EXPROFILEDB } from '../../../../../modules/nf-core/mmseqs/tsv2exprofiledb/main.nf'

workflow test_mmseqs_tsv2exprofiledb {

    db = file('s3://nf-proteinfold-tests/modules_mmseqs_test_data/mmseqs.tar.gz', checkIfExists: true)

    ch_db = UNTAR ( [[:], db] ).untar
        .map { it[1] }
    MMSEQS_TSV2EXPROFILEDB ( ch_db, '' )
}
