#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_ARC_MKGTF } from '../../../../../modules/nf-core/cellrangerarc/mkgtf/main.nf'

workflow test_cellranger_arc_mkgtf {
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    CELLRANGER_ARC_MKGTF ( gtf )
}
