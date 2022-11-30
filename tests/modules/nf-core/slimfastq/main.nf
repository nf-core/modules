#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SLIMFASTQ } from '../../../../modules/nf-core/slimfastq/main.nf'

workflow test_slimfastq_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    SLIMFASTQ ( input )
}

workflow test_slimfastq_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
         file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)]
    ]

    SLIMFASTQ ( input )
}

workflow test_slimfastq_nanopore {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)
    ]

    SLIMFASTQ ( input )
}

workflow test_slimfastq_pacbio {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['ccs_fq_gz'], checkIfExists: true)
    ]

    SLIMFASTQ ( input )
}
