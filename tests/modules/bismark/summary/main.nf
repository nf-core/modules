#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION    } from '../../../../modules/bismark/genomepreparation/main.nf'
include { BISMARK_ALIGN                } from '../../../../modules/bismark/align/main.nf'
include { BISMARK_DEDUPLICATE          } from '../../../../modules/bismark/deduplicate/main.nf'
include { BISMARK_METHYLATIONEXTRACTOR } from '../../../../modules/bismark/methylationextractor/main.nf'
include { BISMARK_SUMMARY              } from '../../../../modules/bismark/summary/main.nf'

workflow test_bismark_summary {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_methylated_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_methylated_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISMARK_GENOMEPREPARATION ( fasta )
    BISMARK_ALIGN ( input, BISMARK_GENOMEPREPARATION.out.index )
    BISMARK_DEDUPLICATE ( BISMARK_ALIGN.out.bam )
    BISMARK_METHYLATIONEXTRACTOR ( BISMARK_DEDUPLICATE.out.bam, BISMARK_GENOMEPREPARATION.out.index )
    BISMARK_SUMMARY (
        BISMARK_ALIGN.out.bam.collect{ meta, bam -> bam },
        BISMARK_ALIGN.out.report.collect{ meta, report -> report },
        BISMARK_DEDUPLICATE.out.report.collect{ meta, bam -> bam },
        BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ meta, report -> report },
        BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ meta, mbias -> mbias }
    )
}
