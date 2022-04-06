#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_FASTQTOSAM } from '../../../../modules/gatk4/fastqtosam/main.nf'

workflow test_gatk4_fastqtosam_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ]

    GATK4_FASTQTOSAM ( input )
}

workflow test_gatk4_fastqtosam_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
            ]

    GATK4_FASTQTOSAM ( input )
}
