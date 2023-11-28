#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPACERANGER_MKREF } from '../../../../../modules/nf-core/spaceranger/mkref/main.nf'

workflow test_spaceranger_mkref {

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    SPACERANGER_MKREF ( fasta,
                        gtf,
                        reference_name )
}
