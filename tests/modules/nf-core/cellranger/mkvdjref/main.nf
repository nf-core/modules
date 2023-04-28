#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKVDJREF } from '../../../../../modules/nf-core/cellranger/mkvdjref/main.nf'

workflow test_cellranger_mkvdjref {

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    CELLRANGER_MKVDJREF ( fasta,
                        gtf,
                        reference_name )
}
