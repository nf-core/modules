#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_INDEX } from '../../../../software/salmon/index/main.nf' addParams( options: [:] )
include { SALMON_QUANT } from '../../../../software/salmon/quant/main.nf' addParams( options: [args: '--minAssignedFrags 1'] )

workflow test_salmon_quant_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true) 
            ]
    genome        = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    transcriptome = file("${launchDir}/tests/data/genomics/sarscov2/genome/transcriptome.fasta", checkIfExists: true)
    gtf           = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.gtf", checkIfExists: true)

    SALMON_INDEX ( genome, transcriptome )
    SALMON_QUANT ( input, SALMON_INDEX.out.index, gtf, transcriptome, false )
}

workflow test_salmon_quant_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_2.fastq.gz", checkIfExists: true) ] 
            ]
    genome        = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    transcriptome = file("${launchDir}/tests/data/genomics/sarscov2/genome/transcriptome.fasta", checkIfExists: true)
    gtf           = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.gtf", checkIfExists: true)

    SALMON_INDEX ( genome, transcriptome )
    SALMON_QUANT ( input, SALMON_INDEX.out.index, gtf, transcriptome, false )
}
