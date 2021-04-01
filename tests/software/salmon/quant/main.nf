#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_INDEX } from '../../../../software/salmon/index/main.nf' addParams( options: [:] )
include { SALMON_QUANT } from '../../../../software/salmon/quant/main.nf' addParams( options: [args: '--minAssignedFrags 1'] )

workflow test_salmon_quant_single_end {

    input = [ [ id:'test', single_end:true ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) 
            ]
    genome_fasta     = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    transcript_fasta = file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)
    gtf              = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    SALMON_INDEX ( genome_fasta, transcript_fasta )
    SALMON_QUANT ( input, SALMON_INDEX.out.index, gtf, transcript_fasta, false )

}

workflow test_salmon_quant_paired_end {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    genome_fasta     = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    transcript_fasta = file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)
    gtf              = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    SALMON_INDEX ( genome_fasta, transcript_fasta )
    SALMON_QUANT ( input, SALMON_INDEX.out.index, gtf, transcript_fasta, false )

}
