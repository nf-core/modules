#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HISAT2_EXTRACTSPLICESITES } from '../../../../software/hisat2/extractsplicesites/main.nf' addParams( options: [:] )
include { HISAT2_BUILD } from '../../../../software/hisat2/build/main.nf' addParams( options: [:] )

workflow test_hisat2_build {
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    HISAT2_EXTRACTSPLICESITES ( gtf )
    HISAT2_BUILD ( fasta, gtf, HISAT2_EXTRACTSPLICESITES.out.txt )
}
