#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION    } from '../../../../software/bismark/genomepreparation/main.nf'    addParams( options: [:] )
include { BISMARK_ALIGN                } from '../../../../software/bismark/align/main.nf'                addParams( options: [:] )
include { BISMARK_DEDUPLICATE          } from '../../../../software/bismark/deduplicate/main.nf'          addParams( options: [:] )
include { BISMARK_METHYLATIONEXTRACTOR } from '../../../../software/bismark/methylationextractor/main.nf' addParams( options: [:] )
include { BISMARK_SUMMARY              } from '../../../../software/bismark/summary/main.nf'              addParams( options: [:] )

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
