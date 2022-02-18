#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUALIMAP_BAMQC } from '../../../../modules/qualimap/bamqc/main.nf'

workflow test_qualimap_bamqc {
    input   = [ [ id:'test', single_end:false ], // meta map
                [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
              ]
    gff     = []

    QUALIMAP_BAMQC ( input, gff )
}
