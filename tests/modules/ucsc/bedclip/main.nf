#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_BEDCLIP } from '../../../../modules/ucsc/bedclip/main.nf'

workflow test_ucsc_bedclip {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_bedgraph'], checkIfExists: true)
            ]
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    UCSC_BEDCLIP ( input, sizes )
}
