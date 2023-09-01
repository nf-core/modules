#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_DEPTH } from '../../../../../modules/nf-core/samtools/depth/main.nf'

workflow test_samtools_depth {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true) ]

    intervals = [ [ id:'bed' ],
                  file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]

    SAMTOOLS_DEPTH ( input, intervals )
}
