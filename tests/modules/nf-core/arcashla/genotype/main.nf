#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARCASHLA_GENOTYPE } from '../../../../../modules/nf-core/arcashla/genotype/main.nf'

workflow test_arcashla_genotype {

    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test2_1_fastq_gz'], checkIfExists: true) ]
    ]

    ARCASHLA_GENOTYPE ( input )
}
