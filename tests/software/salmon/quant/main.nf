#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_INDEX } from '../../../../software/salmon/index/main.nf' addParams( options: [:] )
include { SALMON_QUANT } from '../../../../software/salmon/quant/main.nf' addParams( options: [args: '--minAssignedFrags 1'] )

workflow test_salmon_quant_single_end {

    def genome_fasta     = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true)
    def transcript_fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_cds_from_genomic.fasta", checkIfExists: true)
    def gtf   = file("${launchDir}/tests/data/genomics/sarscov2/gtf/test_genomic.gtf", checkIfExists: true)
    def input = [ [ id:'test', single_end:true ], // meta map
                  file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true) ]

    SALMON_INDEX ( genome_fasta, transcript_fasta )
    SALMON_QUANT ( input, SALMON_INDEX.out.index, gtf, transcript_fasta, false )

}

workflow test_salmon_quant_paired_end {

    def genome_fasta     = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true)
    def transcript_fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_cds_from_genomic.fasta", checkIfExists: true)
    def gtf   = file("${launchDir}/tests/data/genomics/sarscov2/gtf/test_genomic.gtf", checkIfExists: true)
    def input = [ [ id:'test', single_end:false ], // meta map
                  [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true),
                    file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_2.fastq.gz", checkIfExists: true) ] ]

    SALMON_INDEX ( genome_fasta, transcript_fasta )
    SALMON_QUANT ( input, SALMON_INDEX.out.index, gtf, transcript_fasta, false )

}
