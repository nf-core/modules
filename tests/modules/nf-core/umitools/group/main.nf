#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UMITOOLS_EXTRACT                               } from '../../../../../modules/nf-core/umitools/extract/main.nf'
include { BWA_INDEX                                      } from '../../../../../modules/nf-core/bwa/index/main.nf'
include { BWA_MEM                                        } from '../../../../../modules/nf-core/bwa/mem/main.nf'
include { SAMTOOLS_INDEX                                 } from '../../../../../modules/nf-core/samtools/index/main.nf'
include { UMITOOLS_GROUP         } from '../../../../../modules/nf-core/umitools/group/main.nf'

//
// Test with no UMI
//
workflow test_umitools_group_no_umi {
    input = [
        [ id:'test'], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    create_bam     = true
    get_group_info = true

    UMITOOLS_GROUP ( input, create_bam, get_group_info )
}

//
// Test with single-end data with BAM false and group info true
//
workflow test_umitools_group_single_end_info {
    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    fasta = [
        [ id:'sarscov2'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    create_bam     = false
    get_group_info = true

    UMITOOLS_EXTRACT ( input )
    BWA_INDEX ( fasta )
    BWA_MEM ( UMITOOLS_EXTRACT.out.reads, BWA_INDEX.out.index, true )
    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    UMITOOLS_GROUP ( BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]), create_bam, get_group_info )
}

//
// Test with paired-end data with BAM true and group info false
//
workflow test_umitools_group_paired_end_bam {
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
    create_bam     = true
    get_group_info = false

    UMITOOLS_EXTRACT ( input )
    BWA_INDEX ( fasta )
    BWA_MEM ( UMITOOLS_EXTRACT.out.reads, BWA_INDEX.out.index, true )
    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    UMITOOLS_GROUP ( BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]), create_bam, get_group_info )
}

//
// Test with paired-end data BAM true and group info true
//
workflow test_umitools_group_paired_bam_info {
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
    create_bam     = true
    get_group_info = true

    UMITOOLS_EXTRACT ( input )
    BWA_INDEX ( fasta )
    BWA_MEM ( UMITOOLS_EXTRACT.out.reads, BWA_INDEX.out.index, true )
    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    UMITOOLS_GROUP ( BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]), create_bam, get_group_info )
}
