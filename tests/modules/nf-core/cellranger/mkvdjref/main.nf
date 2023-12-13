#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKVDJREF } from '../../../../../modules/nf-core/cellranger/mkvdjref/main.nf'

workflow test_cellranger_mkvdjref {

    reference_name = "homo_sapiens_chr22_reference"

    CELLRANGER_MKVDJREF ( reference_name )
}
