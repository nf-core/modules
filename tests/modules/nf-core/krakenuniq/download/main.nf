#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRAKENUNIQ_DOWNLOAD } from '../../../../../modules/nf-core/krakenuniq/download/main.nf'

workflow test_krakenuniq_download {

    input = "taxonomy"

    KRAKENUNIQ_DOWNLOAD ( input )
}
