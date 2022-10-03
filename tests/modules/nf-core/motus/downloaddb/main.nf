#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOTUS_DOWNLOADDB } from '../../../../../modules/nf-core/motus/downloaddb/main.nf'

workflow test_motus_downloaddb {

    input = file('https://raw.githubusercontent.com/motu-tool/mOTUs/master/motus/downloadDB.py')

    MOTUS_DOWNLOADDB ( input )
}
