#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { UCSC_BEDGRAPHTOBIGWIG  } from "$moduleDir/modules/nf-core/ucsc/bedgraphtobigwig/main.nf"

workflow test_ucsc_bedgraphtobigwig {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_bedgraph'], checkIfExists: true) ]
            ]
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    UCSC_BEDGRAPHTOBIGWIG ( input, sizes )
}
