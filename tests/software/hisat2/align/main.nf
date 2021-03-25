#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HISAT2_EXTRACTSPLICESITES } from '../../../../software/hisat2/extractsplicesites/main.nf' addParams( options: [:] )
include { HISAT2_BUILD } from '../../../../software/hisat2/build/main.nf' addParams( options: [:] )
include { HISAT2_ALIGN } from '../../../../software/hisat2/extractsplicesites/main.nf' addParams( options: [:] )

workflow test_hisat2_align_se {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    HISAT2_EXTRACTSPLICESITES ( gtf )
    HISAT2_BUILD ( fasta, gtf, HISAT2_EXTRACTSPLICESITES.out.txt )
    HISAT2_ALIGN ( input, HISAT2_BUILD, HISAT2_EXTRACTSPLICESITES.out.txt )
}

workflow test_hisat2_align_pe {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
              [ file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    HISAT2_EXTRACTSPLICESITES ( gtf )
    HISAT2_BUILD ( fasta, gtf, HISAT2_EXTRACTSPLICESITES.out.txt )
    HISAT2_ALIGN ( input, HISAT2_BUILD, HISAT2_EXTRACTSPLICESITES.out.txt )
}
