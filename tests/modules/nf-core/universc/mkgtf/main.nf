#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNIVERSC_MKGTF } from '../../../../../modules/nf-core/universc/mkgtf/main.nf'

workflow test_universc_mkgtf {
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    UNIVERSC_MKGTF ( gtf )
}
