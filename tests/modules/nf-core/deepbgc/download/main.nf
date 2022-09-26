#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPBGC_DOWNLOAD } from '../../../../modules/deepbgc/download/main.nf'

workflow test_deepbgc_download {

    DEEPBGC_DOWNLOAD (  )
}
