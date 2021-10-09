#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRATOOLS_PREFETCH } from '../../../../modules/sratools/prefetch/main.nf' addParams( options: [:] )

workflow test_sratools_prefetch {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    SRATOOLS_PREFETCH ( input )
}
