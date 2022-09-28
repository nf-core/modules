#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FASTQ } from '../../../../modules/samtools/fastq/main.nf'

workflow test_samtools_fastq_paired_end_bam {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]

    SAMTOOLS_FASTQ ( input, [], false )
}

workflow test_samtools_fastq_paired_end_cram {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
            ]

    SAMTOOLS_FASTQ ( input,file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true), false )
}

workflow test_samtools_fastq_interleaved {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
            ]

    SAMTOOLS_FASTQ ( input, [], true )
}

workflow test_samtools_fastq_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
            ]

    SAMTOOLS_FASTQ ( input, [], false )
}
