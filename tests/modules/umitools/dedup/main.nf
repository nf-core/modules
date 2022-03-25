#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UMITOOLS_EXTRACT } from '../../../../modules/umitools/extract/main.nf'
include { BWA_INDEX } from '../../../../modules/bwa/index/main.nf'
include { BWA_MEM   } from '../../../../modules/bwa/mem/main.nf'
include { SAMTOOLS_INDEX   } from '../../../../modules/samtools/index/main.nf'
include { UMITOOLS_DEDUP } from '../../../../modules/umitools/dedup/main.nf'

//
// Test with no UMI
//
workflow test_umitools_dedup_no_umi {
    input = [ [ id:'test'], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
            ]

    UMITOOLS_DEDUP ( input )
}

//
// Test with single-end data
//
workflow test_umitools_dedup_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    UMITOOLS_EXTRACT ( input )
    BWA_INDEX ( fasta )
    BWA_MEM ( UMITOOLS_EXTRACT.out.reads, BWA_INDEX.out.index, true )
    SAMTOOLS_INDEX (BWA_MEM.out.bam)
    UMITOOLS_DEDUP(BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]))
}

//
// Test with paired-end data
//
workflow test_umitools_dedup_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    UMITOOLS_EXTRACT ( input )
    BWA_INDEX ( fasta )
    BWA_MEM ( UMITOOLS_EXTRACT.out.reads, BWA_INDEX.out.index, true )
    SAMTOOLS_INDEX (BWA_MEM.out.bam)
    UMITOOLS_DEDUP(BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]))
}
