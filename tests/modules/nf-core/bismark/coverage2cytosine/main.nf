#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION    } from '../../../../../modules/nf-core/bismark/genomepreparation/main.nf'
include { BISMARK_ALIGN                } from '../../../../../modules/nf-core/bismark/align/main.nf'
include { BISMARK_METHYLATIONEXTRACTOR } from '../../../../../modules/nf-core/bismark/methylationextractor/main.nf'
include { BISMARK_COVERAGE2CYTOSINE }    from '../../../../../modules/nf-core/bismark/coverage2cytosine/main.nf'

workflow test_bismark_coverage2cytosine {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_methylated_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_methylated_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISMARK_GENOMEPREPARATION ( fasta )
    BISMARK_ALIGN ( input, BISMARK_GENOMEPREPARATION.out.index )
    BISMARK_METHYLATIONEXTRACTOR ( BISMARK_ALIGN.out.bam, BISMARK_GENOMEPREPARATION.out.index )
    BISMARK_COVERAGE2CYTOSINE ( BISMARK_METHYLATIONEXTRACTOR.out.coverage, BISMARK_GENOMEPREPARATION.out.index )
}
