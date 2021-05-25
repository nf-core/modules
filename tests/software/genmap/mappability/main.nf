#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMAP_INDEX       } from '../../../../software/genmap/index/main.nf'       addParams( options: [:] )
include { GENMAP_MAPPABILITY } from '../../../../software/genmap/mappability/main.nf' addParams( options: [args : '-K 50 -E 2 -w -t -bg'] )

workflow test_genmap_map {

    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GENMAP_INDEX ( input )
    GENMAP_MAPPABILITY ( GENMAP_INDEX.out.index )
}
