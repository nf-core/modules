#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUALIMAP_RNASEQ } from '../../../../../modules/nf-core/qualimap/rnaseq/main.nf'

workflow test_qualimap_rnaseq {
    input   = [ [ id:'test', single_end:false ], // meta map
                [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
              ]
    gff     = []

    QUALIMAP_BAMQC ( input, gff )
}
