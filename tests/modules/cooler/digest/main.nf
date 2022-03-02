#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_DIGEST } from '../../../../modules/cooler/digest/main.nf'

workflow test_cooler_digest {

    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)
    enzyme = "CviQI"

    COOLER_DIGEST ( input, sizes, enzyme )
}
