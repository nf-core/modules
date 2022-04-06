#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKGTF } from '../../../../modules/cellranger/mkgtf/main.nf'

workflow test_cellranger_mkgtf {
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    CELLRANGER_MKGTF ( gtf )
}
