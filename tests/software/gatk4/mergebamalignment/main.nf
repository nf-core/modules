#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MERGEBAMALIGNMENT } from '../../../../software/gatk4/mergebamalignment/main.nf' addParams( options: [:] )

workflow test_gatk4_mergebamalignment {
    input    = [ [ id:'test' ], // meta map
                 file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_single_end.bam", checkIfExists: true) 
               ]
    unmapped = file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_unaligned.bam", checkIfExists: true)
    fasta    = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    dict     = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.dict", checkIfExists: true)

    GATK4_MERGEBAMALIGNMENT ( input, unmapped, fasta, dict )
}
