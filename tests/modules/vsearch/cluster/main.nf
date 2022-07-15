#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VSEARCH_CLUSTER as VSEARCH_CLUSTER_FAST } from '../../../../modules/vsearch/cluster/main.nf'
include { VSEARCH_CLUSTER as VSEARCH_CLUSTER_SIZE } from '../../../../modules/vsearch/cluster/main.nf'
include { VSEARCH_CLUSTER as VSEARCH_CLUSTER_SMALLMEM } from '../../../../modules/vsearch/cluster/main.nf'
include { VSEARCH_CLUSTER as VSEARCH_CLUSTER_UNOISE } from '../../../../modules/vsearch/cluster/main.nf'
include { VSEARCH_CLUSTER as VSEARCH_CLUSTER_USEROUT } from '../../../../modules/vsearch/cluster/main.nf'

workflow test_vsearch_cluster_fast {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    // Test default parameters

    VSEARCH_CLUSTER_FAST ( input )

}

workflow test_vsearch_cluster_size {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    // Test also sam to bam conversion with ext.arg.outoption '--samout'

    VSEARCH_CLUSTER_SIZE ( input )

}

workflow test_vsearch_cluster_smallmem {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    VSEARCH_CLUSTER_SMALLMEM ( input )

}

workflow test_vsearch_cluster_unoise {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    VSEARCH_CLUSTER_UNOISE ( input)

}

workflow test_vsearch_cluster_userout {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    // Test ext.args.outoption "--userout" and ext.args.user_columns

    VSEARCH_CLUSTER_USEROUT ( input )
}
