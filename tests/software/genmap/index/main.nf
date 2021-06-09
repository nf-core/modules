#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMAP_INDEX } from '../../../../software/genmap/index/main.nf' addParams( options: [:] )

workflow test_genmap_index {

    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GENMAP_INDEX ( input )
}
