#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPBGC_DOWNLOAD } from '../../../../modules/deepbgc/download/main.nf'

workflow test_deepbgc_download {

    //input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    DEEPBGC_DOWNLOAD (  )
}
