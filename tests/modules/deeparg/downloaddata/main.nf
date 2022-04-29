#!/usr/bin/env nextflow



include { DEEPARG_DOWNLOADDATA } from '../../../../modules/deeparg/downloaddata/main.nf'

workflow test_deeparg_downloaddata {
    DEEPARG_DOWNLOADDATA ( )
}
