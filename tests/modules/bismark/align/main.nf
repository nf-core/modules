#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION         } from '../../../../modules/bismark/genomepreparation/main.nf'
include { BISMARK_ALIGN as BISMARK_ALIGN_SE } from '../../../../modules/bismark/align/main.nf'
include { BISMARK_ALIGN as BISMARK_ALIGN_PE } from '../../../../modules/bismark/align/main.nf'

//
// Test with single-end data
//
workflow test_bismark_align_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_methylated_1_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISMARK_GENOMEPREPARATION ( fasta )
    BISMARK_ALIGN_SE ( input, BISMARK_GENOMEPREPARATION.out.index )
}

//
// Test with paired-end data
//
workflow test_bismark_align_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_methylated_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_methylated_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISMARK_GENOMEPREPARATION ( fasta )
    BISMARK_ALIGN_PE ( input, BISMARK_GENOMEPREPARATION.out.index )
}
