#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GETFASTA } from '../../../../modules/bedtools/getfasta/main.nf'

workflow test_bedtools_getfasta {
    bed   = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    BEDTOOLS_GETFASTA ( bed, fasta )
}
