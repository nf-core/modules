#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_INDEX } from '../../../../modules/bbmap/index/main.nf'
include { BBMAP_ALIGN } from '../../../../modules/bbmap/align/main.nf'
include { BBMAP_ALIGN as BBMAP_ALIGN_PIGZ } from '../../../../modules/bbmap/align/main.nf'

workflow test_bbmap_align_paired_end_fasta_ref {

    input = [ [ id:'test', single_end:false ], // meta map
                [
                    file( params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file( params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BBMAP_ALIGN ( input, fasta )
}

workflow test_bbmap_align_paired_end_index_ref {

    input = [ [ id:'test', single_end:false ], // meta map
                [
                    file( params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file( params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BBMAP_INDEX ( fasta )
    BBMAP_ALIGN ( input, BBMAP_INDEX.out.index )
}

workflow test_bbmap_align_single_end_index_ref {

    input = [ [ id:'test', single_end:true ], // meta map
                file( params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BBMAP_INDEX ( fasta )
    BBMAP_ALIGN ( input, BBMAP_INDEX.out.index )
}

workflow test_bbmap_align_paired_end_index_ref_pigz {

    input = [ [ id:'test', single_end:false ], // meta map
                [
                    file( params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file( params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BBMAP_INDEX ( fasta )
    BBMAP_ALIGN_PIGZ ( input, BBMAP_INDEX.out.index )
}
