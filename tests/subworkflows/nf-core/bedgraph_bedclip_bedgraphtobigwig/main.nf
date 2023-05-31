#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG } from '../../../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'

workflow test_bedgraph_bedclip_bedgraphtobigwig {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_bedgraph'], checkIfExists: true)
            ]
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG ( input, sizes )
}
