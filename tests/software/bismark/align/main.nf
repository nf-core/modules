#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION         } from '../../../../software/bismark/genomepreparation/main.nf' addParams( options: [:]                               )
include { BISMARK_ALIGN as BISMARK_ALIGN_SE } from '../../../../software/bismark/align/main.nf'             addParams( options: [ publish_dir:'test_single_end' ] )
include { BISMARK_ALIGN as BISMARK_ALIGN_PE } from '../../../../software/bismark/align/main.nf'             addParams( options: [ publish_dir:'test_paired_end' ] )

/*
 * Test with single-end data
 */
workflow test_bismark_align_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_methylated_1.fastq.gz", checkIfExists: true) ]
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    BISMARK_GENOMEPREPARATION ( fasta )
    BISMARK_ALIGN_SE ( input, BISMARK_GENOMEPREPARATION.out.index )
}

/*
 * Test with paired-end data
 */
workflow test_bismark_align_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_methylated_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_methylated_2.fastq.gz", checkIfExists: true) ]
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    BISMARK_GENOMEPREPARATION ( fasta )
    BISMARK_ALIGN_PE ( input, BISMARK_GENOMEPREPARATION.out.index )
}
