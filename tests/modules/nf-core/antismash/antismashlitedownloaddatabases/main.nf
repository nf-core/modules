#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR as UNTAR1 } from '../../../../../modules/nf-core/untar/main.nf'
include { UNTAR as UNTAR2 } from '../../../../../modules/nf-core/untar/main.nf'
include { UNTAR as UNTAR3 } from '../../../../../modules/nf-core/untar/main.nf'
include { ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES } from '../../../../../modules/nf-core/antismash/antismashlitedownloaddatabases/main.nf'

workflow test_antismash_antismashlitedownloaddatabases {

    // hash 91bb8781c576967e23d2c5315dd4d43213575033 is for antismash version 6.1.1
    antismash_css = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/91bb8781c576967e23d2c5315dd4d43213575033/data/delete_me/antismash/css.tar.gz', checkIfExists: true)
    ]

    antismash_detection = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/91bb8781c576967e23d2c5315dd4d43213575033/data/delete_me/antismash/detection.tar.gz', checkIfExists: true)
    ]

    antismash_modules = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/91bb8781c576967e23d2c5315dd4d43213575033/data/delete_me/antismash/modules.tar.gz', checkIfExists: true)
    ]


    UNTAR1 ( antismash_css )
    UNTAR2 ( antismash_detection )
    UNTAR3 ( antismash_modules )

    ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES ( UNTAR1.out.untar.map{ it[1] },  UNTAR2.out.untar.map{ it[1] }, UNTAR3.out.untar.map{ it[1] } )
}
