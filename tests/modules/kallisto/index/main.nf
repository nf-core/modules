#!/usr/bin/env nextflow



include { KALLISTO_INDEX } from '../../../../modules/kallisto/index/main.nf'

workflow test_kallisto_index {

    def input = []
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    KALLISTO_INDEX ( input )
}
