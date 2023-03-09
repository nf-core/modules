#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UMITOOLS_EXTRACT } from '../../../../../modules/nf-core/umitools/extract/main.nf'
include { BWA_INDEX        } from '../../../../../modules/nf-core/bwa/index/main.nf'
include { BWA_MEM          } from '../../../../../modules/nf-core/bwa/mem/main.nf'
include { SAMTOOLS_INDEX   } from '../../../../../modules/nf-core/samtools/index/main.nf'
include { UMITOOLS_DEDUP   } from '../../../../../modules/nf-core/umitools/dedup/main.nf'

//
// Test with no UMI
//
workflow test_umitools_dedup_no_umi {
    input = [ 
        [ id:'test'], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    get_output_stats = false

    UMITOOLS_DEDUP ( input, get_output_stats )
}

//
// Test with single-end data without --output-stats
//
workflow test_umitools_dedup_single_end_no_stats {
    input = [ 
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    fasta = [
        [ id:'sarscov2'], 
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    get_output_stats = false

    UMITOOLS_EXTRACT ( input )
    BWA_INDEX ( fasta )
    BWA_MEM ( UMITOOLS_EXTRACT.out.reads, BWA_INDEX.out.index, true )
    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    UMITOOLS_DEDUP ( BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]), get_output_stats )
}

//
// Test with paired-end data without --output-stats
//
workflow test_umitools_dedup_paired_end_no_stats {
    input = [ 
        [ id:'test', single_end:false ], // meta map
        [ 
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) 
        ]
    ]
    fasta = [
        [ id:'sarscov2'], 
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    get_output_stats = false

    UMITOOLS_EXTRACT ( input )
    BWA_INDEX ( fasta )
    BWA_MEM ( UMITOOLS_EXTRACT.out.reads, BWA_INDEX.out.index, true )
    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    UMITOOLS_DEDUP ( BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]), get_output_stats )
}

//
// Test with paired-end data with --output-stats
//
workflow test_umitools_dedup_paired_end_stats {
    input = [ 
        [ id:'test', single_end:false ], // meta map
        [ 
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) 
        ]
    ]
    fasta = [
        [ id:'sarscov2'], 
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    get_output_stats = true

    UMITOOLS_EXTRACT ( input )
    BWA_INDEX ( fasta )
    BWA_MEM ( UMITOOLS_EXTRACT.out.reads, BWA_INDEX.out.index, true )
    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    UMITOOLS_DEDUP ( BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]), get_output_stats )
}
