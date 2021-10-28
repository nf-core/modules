#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CIRCULARMAPPER_CIRCULARMAPPER } from '../../../../modules/circularmapper/circularmapper/main.nf' addParams( options: [:] )

workflow test_circularmapper_circularmapper {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    CIRCULARMAPPER_CIRCULARMAPPER ( input )
}
