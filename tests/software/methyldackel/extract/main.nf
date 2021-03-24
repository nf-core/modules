#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METHYLDACKEL_EXTRACT } from '../../../../software/methyldackel/extract/main.nf' addParams( options: [:] )

workflow test_methyldackel_extract {
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_methylated_paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_methylated_paired_end.sorted.bam.bai", checkIfExists: true) ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    fai   = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta.fai", checkIfExists: true)

    METHYLDACKEL_EXTRACT ( input, fasta, fai )
}
