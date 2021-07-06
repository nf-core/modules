#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFFREAD } from '../../../software/gffread/main.nf' addParams( options: [:] )

workflow test_gffread {
    input = file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true)

    GFFREAD ( input )
}
