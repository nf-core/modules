#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../modules/bwa/index/main.nf' addParams( options: [:] )
include { BWA_MEM } from '../../../../modules/bwa/mem/main.nf' addParams( options: [:] )

//
// Test with single-end data
//
workflow test_bwa_mem_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWA_INDEX ( fasta )
    BWA_MEM ( input, BWA_INDEX.out.index )
}

//
// Test with paired-end data
//
workflow test_bwa_mem_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWA_INDEX ( fasta )
    BWA_MEM ( input, BWA_INDEX.out.index )
}
