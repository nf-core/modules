#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNIVERSC_CELLRANGER_OS_MKGTF } from '../../../../../modules/nf-core/universc/cellranger_os_mkgtf/main.nf'

workflow test_universc_mkgtf {
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    UNIVERSC_CELLRANGER_OS_MKGTF ( gtf )
}
