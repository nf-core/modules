#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION    } from '../../../../modules/bismark/genomepreparation/main.nf'
include { BISMARK_ALIGN                } from '../../../../modules/bismark/align/main.nf'
include { BISMARK_DEDUPLICATE          } from '../../../../modules/bismark/deduplicate/main.nf'
include { BISMARK_METHYLATIONEXTRACTOR } from '../../../../modules/bismark/methylationextractor/main.nf'
include { BISMARK_REPORT               } from '../../../../modules/bismark/report/main.nf'

workflow test_bismark_report {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_methylated_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_methylated_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISMARK_GENOMEPREPARATION ( fasta )
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
