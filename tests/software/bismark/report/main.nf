#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION    } from '../../../../software/bismark/genomepreparation/main.nf' addParams( options: [:] )
include { BISMARK_ALIGN                 } from '../../../../software/bismark/align/main.nf' addParams( options: [:] )
include { BISMARK_DEDUPLICATE           } from '../../../../software/bismark/deduplicate/main.nf' addParams( options: [:] )
include { BISMARK_METHYLATIONEXTRACTOR } from '../../../../software/bismark/methylationextractor/main.nf' addParams( options: [:] )
include { BISMARK_REPORT                } from '../../../../software/bismark/report/main.nf' addParams( options: [:] )

workflow test_bismark_report {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_methylated_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_methylated_2.fastq.gz", checkIfExists: true) ] ]

    BISMARK_GENOMEPREPARATION ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )
    BISMARK_ALIGN ( input, BISMARK_GENOMEPREPARATION.out.index )
    BISMARK_DEDUPLICATE ( BISMARK_ALIGN.out.bam )
    BISMARK_METHYLATIONEXTRACTOR ( BISMARK_DEDUPLICATE.out.bam, BISMARK_GENOMEPREPARATION.out.index )

    BISMARK_REPORT (
        BISMARK_ALIGN.out.report
            .join(BISMARK_DEDUPLICATE.out.report)
            .join(BISMARK_METHYLATIONEXTRACTOR.out.report)
            .join(BISMARK_METHYLATIONEXTRACTOR.out.mbias)
    )
}
