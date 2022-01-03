#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION    } from '../../../../modules/bismark/genomepreparation/main.nf'
include { BISMARK_METHYLATIONEXTRACTOR } from '../../../../modules/bismark/methylationextractor/main.nf'

workflow test_bismark_methylationextractor {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_bam'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISMARK_GENOMEPREPARATION ( fasta )
    BISMARK_METHYLATIONEXTRACTOR ( input, BISMARK_GENOMEPREPARATION.out.index )
}
