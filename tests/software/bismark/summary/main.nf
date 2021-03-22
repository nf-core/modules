#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOME_PREPARATION    } from '../../../../software/bismark/genome_preparation/main.nf' addParams( options: [:] )
include { BISMARK_ALIGN                 } from '../../../../software/bismark/align/main.nf' addParams( options: [:] )
include { BISMARK_DEDUPLICATE           } from '../../../../software/bismark/deduplicate/main.nf' addParams( options: [:] )
include { BISMARK_METHYLATION_EXTRACTOR } from '../../../../software/bismark/methylation_extractor/main.nf' addParams( options: [:] )
include { BISMARK_SUMMARY               } from '../../../../software/bismark/summary/main.nf' addParams( options: [:] )

workflow test_bismark_summary {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_methylated_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_methylated_2.fastq.gz", checkIfExists: true) ] ]

    BISMARK_GENOME_PREPARATION ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )
    BISMARK_ALIGN ( input, BISMARK_GENOME_PREPARATION.out.index )
    BISMARK_DEDUPLICATE ( BISMARK_ALIGN.out.bam )
    BISMARK_METHYLATION_EXTRACTOR ( BISMARK_DEDUPLICATE.out.bam, BISMARK_GENOME_PREPARATION.out.index )

    BISMARK_SUMMARY (
        BISMARK_ALIGN.out.bam.collect{ meta, bam -> bam },
        BISMARK_ALIGN.out.report.collect{ meta, report -> report },
        BISMARK_DEDUPLICATE.out.bam.collect{ meta, bam -> bam },
        BISMARK_METHYLATION_EXTRACTOR.out.report.collect{ meta, report -> report },
        BISMARK_METHYLATION_EXTRACTOR.out.mbias.collect{ meta, mbias -> mbias }
    )
}
