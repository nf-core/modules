#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPARG_DOWNLOADDATA } from '../../../../../modules/nf-core/deeparg/downloaddata/main.nf'

workflow test_deeparg_downloaddata {
    DEEPARG_DOWNLOADDATA ( )
}
