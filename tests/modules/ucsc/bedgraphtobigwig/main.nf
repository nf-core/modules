#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_BEDGRAPHTOBIGWIG  } from '../../../../modules/ucsc/bedgraphtobigwig/main.nf'

workflow test_ucsc_bedgraphtobigwig {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_bedgraph'], checkIfExists: true) ]
            ]
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    UCSC_BEDGRAPHTOBIGWIG ( input, sizes )
}
