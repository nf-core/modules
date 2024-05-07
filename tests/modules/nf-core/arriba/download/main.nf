#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARRIBA_DOWNLOAD } from '../../../../../modules/nf-core/arriba/download/main.nf'

workflow test_arriba_download {

    ARRIBA_DOWNLOAD ()
}
