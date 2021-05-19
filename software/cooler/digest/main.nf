#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_DIGEST } from '../../../../software/cooler/digest/main.nf' addParams( options: [:] )

workflow test_cooler_digest {

    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)
    enzyme = "CviQI"

    COOLER_DIGEST ( input, sizes, enzyme )
}
