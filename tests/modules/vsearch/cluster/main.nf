#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VSEARCH_CLUSTER } from '../../../../modules/vsearch/cluster/main.nf'

workflow test_vsearch_cluster {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    VSEARCH_CLUSTER ( input )
}
