#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HISAT2_EXTRACTSPLICESITES } from '../../../../modules/nf-core/hisat2/extractsplicesites/main.nf'
include { HISAT2_BUILD              } from '../../../../modules/nf-core/hisat2/build/main.nf'
include { FASTQ_ALIGN_HISAT2        } from '../../../../subworkflows/nf-core/fastq_align_hisat2/main.nf'

workflow test_fastq_align_hisat2_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [ id:'test' ],
        [
            file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        ]
    ]
    gtf   = [ [id:'test'],
            file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)
    ]

    HISAT2_EXTRACTSPLICESITES ( gtf )
    HISAT2_BUILD ( fasta, gtf, HISAT2_EXTRACTSPLICESITES.out.txt )
    FASTQ_ALIGN_HISAT2 ( input, HISAT2_BUILD.out.index, HISAT2_EXTRACTSPLICESITES.out.txt, fasta )
}

workflow test_fastq_align_hisat2_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [ id:'test' ],
        [
            file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        ]
    ]
    gtf   = [ [id:'test'],
            file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)
    ]

    HISAT2_EXTRACTSPLICESITES ( gtf )
    HISAT2_BUILD ( fasta, gtf, HISAT2_EXTRACTSPLICESITES.out.txt )
    FASTQ_ALIGN_HISAT2 ( input, HISAT2_BUILD.out.index, HISAT2_EXTRACTSPLICESITES.out.txt, fasta )
}
