#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGERARC_MKGTF } from '../../../../../modules/nf-core/cellrangerarc/mkgtf/main.nf'

workflow test_cellrangerarc_mkgtf {
    gtf = file(params.test_data['mus_musculus']['genome']['genome_19_gtf'], checkIfExists: true)

    CELLRANGERARC_MKGTF ( gtf )
}
