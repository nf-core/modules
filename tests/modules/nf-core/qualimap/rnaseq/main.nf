#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUALIMAP_RNASEQ } from '../../../../../modules/nf-core/qualimap/rnaseq/main.nf'

workflow test_qualimap_rnaseq {
    bam   = [ [ id:'test', single_end:false ], // meta map
                [ file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]
              ]
    gtf = [
        [ id:'test_gtf' ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true) ]
    ]

    QUALIMAP_RNASEQ ( bam, gtf )
}
