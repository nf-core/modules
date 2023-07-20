#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HISAT2_EXTRACTSPLICESITES } from '../../../../../modules/nf-core/hisat2/extractsplicesites/main.nf'
include { HISAT2_BUILD              } from '../../../../../modules/nf-core/hisat2/build/main.nf'
include { HISAT2_ALIGN              } from '../../../../../modules/nf-core/hisat2/align/main.nf'

workflow test_hisat2_align_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [ [id:'genome'],
            file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]
    gtf   = [ [id:'test'],
            file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)
            ]

    HISAT2_EXTRACTSPLICESITES ( gtf )
    HISAT2_BUILD ( fasta, gtf, HISAT2_EXTRACTSPLICESITES.out.txt )
    HISAT2_ALIGN ( input, HISAT2_BUILD.out.index, HISAT2_EXTRACTSPLICESITES.out.txt )
}

workflow test_hisat2_align_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [ [id:'genome'],
            file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]
    gtf   = [ [id:'test'],
            file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)
            ]

    HISAT2_EXTRACTSPLICESITES ( gtf )
    HISAT2_BUILD ( fasta, gtf, HISAT2_EXTRACTSPLICESITES.out.txt )
    HISAT2_ALIGN ( input, HISAT2_BUILD.out.index, HISAT2_EXTRACTSPLICESITES.out.txt )
}

workflow test_hisat2_align_single_end_no_ss {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [ [id:'genome'],
            file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]
    gtf   = [ [id:'test'],
            file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)
            ]

    HISAT2_BUILD ( fasta, [[],[]], [[],[]] )
    HISAT2_ALIGN ( input, [[],[]], [[],[]] )
}