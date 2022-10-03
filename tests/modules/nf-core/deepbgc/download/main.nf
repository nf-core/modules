#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPBGC_DOWNLOAD } from '../../../../../modules/nf-core/deepbgc/download/main.nf'

workflow test_deepbgc_download {

    DEEPBGC_DOWNLOAD (  )
}
