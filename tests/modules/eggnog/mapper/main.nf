#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EGGNOG_MAPPER } from '../../../../modules/eggnog/mapper/main.nf' addParams( options: [:] )

workflow test_eggnog_mapper {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
    db = []

    EGGNOG_MAPPER ( input, db )
}
