#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENOMAD_DOWNLOAD } from '../../../../../modules/nf-core/genomad/download/main.nf'

workflow test_genomad_download {

    GENOMAD_DOWNLOAD ( )
}
