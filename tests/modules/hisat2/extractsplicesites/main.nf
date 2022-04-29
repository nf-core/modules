#!/usr/bin/env nextflow



include { HISAT2_EXTRACTSPLICESITES } from '../../../../modules/hisat2/extractsplicesites/main.nf'

workflow test_hisat2_extractsplicesites {
    gtf = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    HISAT2_EXTRACTSPLICESITES ( gtf )
}
