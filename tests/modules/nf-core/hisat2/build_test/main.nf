#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HISAT2_EXTRACTSPLICESITES } from '../../../../modules/hisat2/extractsplicesites/main.nf'
include { HISAT2_BUILD              } from '../../../../modules/hisat2/build/main.nf'

workflow test_hisat2_build {
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    HISAT2_EXTRACTSPLICESITES ( gtf )
    HISAT2_BUILD ( fasta, gtf, HISAT2_EXTRACTSPLICESITES.out.txt )
}
