#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FAIDX } from '../../../../modules/samtools/faidx/main.nf'

workflow test_samtools_faidx {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    SAMTOOLS_FAIDX ( input )
}
