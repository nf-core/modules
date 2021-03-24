#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ALLELECOUNTER } from '../../../software/allelecounter/main.nf' addParams( options: [:] )

workflow test_allelecounter {
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam.bai", checkIfExists: true)
            ]
    positions = [ file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true) ]

    ALLELECOUNTER ( input, positions )
}
