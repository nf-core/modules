#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GPROFILER2_GOST } from '../../../../../modules/nf-core/gprofiler2/gost/main.nf'

workflow test_gprofiler2_gost {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GPROFILER2_GOST ( input )
}
