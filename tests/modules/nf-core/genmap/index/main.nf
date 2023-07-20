#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMAP_INDEX } from '../../../../../modules/nf-core/genmap/index/main.nf'

workflow test_genmap_index {

    input = [
        [id:"test"],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    GENMAP_INDEX ( input )
}
