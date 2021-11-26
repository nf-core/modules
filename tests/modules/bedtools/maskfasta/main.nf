#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_MASKFASTA } from '../../../../modules/bedtools/maskfasta/main.nf'

workflow test_bedtools_maskfasta {
    bed   = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    BEDTOOLS_MASKFASTA ( bed, fasta )
}
