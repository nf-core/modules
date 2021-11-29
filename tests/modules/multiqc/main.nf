#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC  } from '../../../modules/fastqc/main.nf'
include { MULTIQC } from '../../../modules/multiqc/main.nf'

workflow test_multiqc {
    input = [ [ id: 'test', single_end: false ],
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    FASTQC  ( input )
    MULTIQC ( FASTQC.out.zip.collect { it[1] } )
}
