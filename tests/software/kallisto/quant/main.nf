#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTO_INDEX } from '../../../../software/kallisto/quant/main.nf' addParams( options: [:] )
include { KALLISTO_QUANT } from '../../../../software/kallisto/quant/main.nf' addParams( options: [:] )

workflow test_kallisto_quant_single_end {

    def genome_fasta     = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
    def transcript_fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_transcriptome.fasta", checkIfExists: true)
    def gtf   = file("${launchDir}/tests/data/genomics/sarscov2/gtf/test_genome.gtf", checkIfExists: true)
    def input = [ [ id:'test', single_end:true ], // meta map
                  file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true) ]

    KALLISTO_INDEX ( genome_fasta, transcript_fasta )
    KALLISTO_QUANT ( input, KALLISTO_INDEX.out.index, gtf, transcript_fasta, false )

}

workflow test_kallisto_quant_paired_end {

    def genome_fasta     = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
    def transcript_fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_transcriptome.fasta", checkIfExists: true)
    def gtf   = file("${launchDir}/tests/data/genomics/sarscov2/gtf/test_genome.gtf", checkIfExists: true)
    def input = [ [ id:'test', single_end:false ], // meta map
                  [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true),
                    file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_2.fastq.gz", checkIfExists: true) ] ]

    KALLISTO_INDEX ( genome_fasta, transcript_fasta )
    KALLISTO_QUANT ( input, KALLISTO_INDEX.out.index, gtf, transcript_fasta, false )

}

