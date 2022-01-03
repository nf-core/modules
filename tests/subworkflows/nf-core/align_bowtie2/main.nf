#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../../modules/bowtie2/build/main.nf' addParams( options: [:] )
include { ALIGN_BOWTIE2 } from '../../../../subworkflows/nf-core/align_bowtie2/main.nf' addParams( 'samtools_sort_options': ['suffix': '.sorted'] )

workflow test_align_bowtie2_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BOWTIE2_BUILD ( fasta )
    ALIGN_BOWTIE2 ( input, BOWTIE2_BUILD.out.index )
}

workflow test_align_bowtie2_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BOWTIE2_BUILD ( fasta )
    ALIGN_BOWTIE2 ( input, BOWTIE2_BUILD.out.index )
}
