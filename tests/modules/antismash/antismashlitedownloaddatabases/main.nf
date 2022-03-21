#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES } from '../../../../modules/antismash/antismashlitedownloaddatabases/main.nf'

workflow test_antismash_antismashlitedownloaddatabases {

    input = ''

    ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES ( input )
}
