#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TIDDIT_SV } from '../../../../software/tiddit/sv/main.nf' addParams( options: [:] )

workflow test_tiddit_sv {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) ] 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    fai   = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta.fai", checkIfExists: true)

    TIDDIT_SV ( input, fasta, fai )
}

workflow test_tiddit_sv_no_ref {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) ] 
            ]

    TIDDIT_SV ( input, [], [] )
}
