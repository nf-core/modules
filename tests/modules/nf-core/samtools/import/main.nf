#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_IMPORT } from '../../../../../modules/nf-core/samtools/import/main.nf'

workflow test_samtools_import_single {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    SAMTOOLS_IMPORT ( input )
}

workflow test_samtools_import_paired {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true),
    ]

    SAMTOOLS_IMPORT ( input )
}

workflow test_samtools_import_interleaved {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_interleaved_fastq_gz'], checkIfExists: true)
    ]

    SAMTOOLS_IMPORT ( input )
}
