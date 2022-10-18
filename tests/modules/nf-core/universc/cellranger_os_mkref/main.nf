#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNIVERSC_MKREF } from '../../../../../modules/nf-core/universc/cellranger_os_mkref/main.nf'

workflow test_universc_mkref {

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    UNIVERSC_MKREF ( fasta,
        gtf,
        reference_name )
}
