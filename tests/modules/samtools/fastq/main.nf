#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FASTQ } from '../../../../modules/samtools/fastq/main.nf'

workflow test_samtools_fastq {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]

    SAMTOOLS_FASTQ ( input )
}
