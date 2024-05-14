#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPACERANGER_MKGTF } from '../../../../../modules/nf-core/spaceranger/mkgtf/main.nf'

workflow test_spaceranger_mkgtf {
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    SPACERANGER_MKGTF ( gtf )
}
