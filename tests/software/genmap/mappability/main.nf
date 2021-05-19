#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMAP_MAPPABILITY } from '../../../../software/genmap/mappability/main.nf' addParams( options: [args : '-K 50 -E 2'] )

workflow test_genmap_mappability {

    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GENMAP_MAPPABILITY ( input )
}
