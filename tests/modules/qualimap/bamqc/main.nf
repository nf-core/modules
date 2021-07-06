#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUALIMAP_BAMQC } from '../../../../software/qualimap/bamqc/main.nf' addParams( options: [:] )

workflow test_qualimap_bamqc {
    input   = [ [ id:'test', single_end:false ], // meta map
                [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
              ]
    gff     = file("dummy_file.txt")
    use_gff = false

    QUALIMAP_BAMQC ( input, gff, use_gff )
}
