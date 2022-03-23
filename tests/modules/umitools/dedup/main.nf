#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UMITOOLS_EXTRACT } from '../../../../modules/umitools/extract/main.nf'
include { BWA_INDEX } from '../../../../modules/bwa/index/main.nf'
include { BWA_MEM   } from '../../../../modules/bwa/mem/main.nf'
include { UMITOOLS_DEDUP } from '../../../../modules/umitools/dedup/main.nf'

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
    UMITOOLS_DEDUP(BWA_MEM.out.bam)
}
