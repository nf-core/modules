#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_CUTN } from '../../../../../modules/nf-core/seqtk/cutn/main.nf'

//
// Test with single-end data
//

workflow test_seqtk_cutn {

    input = [ [ id:'test', single_end:true ], // meta map
              file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true) ]

    SEQTK_CUTN ( input )
}