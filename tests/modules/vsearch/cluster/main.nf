#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VSEARCH_CLUSTER } from '../../../../modules/vsearch/cluster/main.nf'
include { VSEARCH_CLUSTER as VSEARCH_CLUSTER_SMALLMEM } from '../../../../modules/vsearch/cluster/main.nf'
include { VSEARCH_CLUSTER as VSEARCH_CLUSTER_UNOISE } from '../../../../modules/vsearch/cluster/main.nf'

workflow test_vsearch_cluster_fast {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    clusteroption = "abcd" // Nonsense text to check default case.
    outoption = "abcd"  // Nonsense text to check default case.

    VSEARCH_CLUSTER ( input, clusteroption, idcutoff, outoption, user_columns )

}

workflow test_vsearch_cluster_size {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    clusteroption = "size"
    outoption = "samout"  // Test also sam to bam conversion

    VSEARCH_CLUSTER ( input, clusteroption, idcutoff, outoption, user_columns )

}

workflow test_vsearch_cluster_smallmem {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    clusteroption = "smallmem"
    outoption = "abcd"  // Nonsense text to check default case.

    VSEARCH_CLUSTER_SMALLMEM ( input, clusteroption, idcutoff, outoption, user_columns )

}

workflow test_vsearch_cluster_unoise {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    clusteroption = "unoise"
    outoption = "abcd"

    VSEARCH_CLUSTER_UNOISE ( input, clusteroption, idcutoff, outoption, user_columns )

}

workflow test_vsearch_cluster_userout {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    clusteroption = "abcd" // Nonsense text to check default case.
    outoption = "userout"

    VSEARCH_CLUSTER ( input, clusteroption, idcutoff, outoption, user_columns )
}
