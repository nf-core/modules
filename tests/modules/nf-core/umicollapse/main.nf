#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UMITOOLS_EXTRACT } from '../../../../modules/nf-core/umitools/extract/main.nf'
include { BWA_INDEX        } from '../../../../modules/nf-core/bwa/index/main.nf'
include { BWA_MEM          } from '../../../../modules/nf-core/bwa/mem/main.nf'
include { SAMTOOLS_INDEX   } from '../../../../modules/nf-core/samtools/index/main.nf'
include { UMICOLLAPSE      } from '../../../../modules/nf-core/umicollapse/main.nf'


//
// Test with single-end data 
//
workflow test_umicollapse_single_end {
    input = [ 
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    fasta = [
        [ id:'sarscov2'], 
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    UMITOOLS_EXTRACT ( input )
    BWA_INDEX ( fasta )
    BWA_MEM ( UMITOOLS_EXTRACT.out.reads, BWA_INDEX.out.index, true )
    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    UMICOLLAPSE ( BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]) )
}

//
// Test with paired-end data 

workflow test_umicollapse_paired_end {
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

    UMITOOLS_EXTRACT ( input )
    BWA_INDEX ( fasta )
    BWA_MEM ( UMITOOLS_EXTRACT.out.reads, BWA_INDEX.out.index, true )
    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    UMICOLLAPSE ( BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]) )
}
