#!/usr/bin/env nextflow



include { GENMAP_INDEX } from '../../../../modules/genmap/index/main.nf'

workflow test_genmap_index {

    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GENMAP_INDEX ( input )
}
