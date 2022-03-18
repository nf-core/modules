#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPARG_DOWNLOADDATA } from '../../../../modules/deeparg/downloaddata/main.nf'

workflow test_deeparg_downloaddata {
    dummy = file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    DEEPARG_DOWNLOADDATA ( dummy )
}
