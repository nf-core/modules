#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FASTQ } from '../../../../software/samtools/fastq/main.nf' addParams( options: [:] )

workflow test_samtools_fastq {
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) 
            ]

    SAMTOOLS_FASTQ ( input )
}
