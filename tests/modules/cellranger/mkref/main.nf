#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKREF } from '../../../../modules/cellranger/mkref/main.nf'

workflow test_cellranger_mkref {

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    CELLRANGER_MKREF ( fasta,
                        gtf,
                        reference_name )
}
