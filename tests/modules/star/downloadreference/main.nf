#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_DOWNLOADREFERENCE } from '../../../../modules/star/downloadreference/main.nf' addParams( options: [:] )

workflow test_star_downloadreference {

    params.genome = ""
    genome = params.genome

    STAR_DOWNLOADREFERENCE ( genome )
}
