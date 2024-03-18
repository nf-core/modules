#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_DICT } from '../../../../../modules/nf-core/samtools/dict/main.nf'

workflow test_samtools_dict {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    SAMTOOLS_DICT ( input )
}
