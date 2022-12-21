#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR as UNTAR1 } from '../../../../../modules/nf-core/untar/main.nf'
include { UNTAR as UNTAR2 } from '../../../../../modules/nf-core/untar/main.nf'
include { UNTAR as UNTAR3 } from '../../../../../modules/nf-core/untar/main.nf'
include { ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES } from '../../../../../modules/nf-core/antismash/antismashlitedownloaddatabases/main.nf'

workflow test_antismash_antismashlitedownloaddatabases {
    // hash 9d816b1327c5f73cd0113af0966197db1956c840 is additional files for for antismash-lite 6.0.1
    input1 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/9d816b1327c5f73cd0113af0966197db1956c840/data/delete_me/antismash/css.tar.gz', checkIfExists: true)
    ]
    input2 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/9d816b1327c5f73cd0113af0966197db1956c840/data/delete_me/antismash/detection.tar.gz', checkIfExists: true)
    ]
    input3 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/9d816b1327c5f73cd0113af0966197db1956c840/data/delete_me/antismash/modules.tar.gz', checkIfExists: true)
    ]

    UNTAR1 ( input1 )
    UNTAR2 ( input2 )
    UNTAR3 ( input3 )

    ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES ( UNTAR1.out.untar.map{ it[1] },  UNTAR2.out.untar.map{ it[1] }, UNTAR3.out.untar.map{ it[1] } )
}
