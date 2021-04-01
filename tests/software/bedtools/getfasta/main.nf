#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GETFASTA } from '../../../../software/bedtools/getfasta/main.nf' addParams( options: [:] )

workflow test_bedtools_getfasta {
    bed   = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    BEDTOOLS_GETFASTA ( bed, fasta )
}
