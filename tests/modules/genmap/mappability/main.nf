#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMAP_INDEX       } from '../../../../modules/genmap/index/main.nf'
include { GENMAP_MAPPABILITY } from '../../../../modules/genmap/mappability/main.nf'

workflow test_genmap_map {

    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GENMAP_INDEX ( input )
    GENMAP_MAPPABILITY ( GENMAP_INDEX.out.index )
}
