#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUALIMAP_BAMQC } from '../../../../software/qualimap/bamqc/main.nf' addParams( options: [:] )

workflow test_qualimap_bamqc {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test-sc2-artic-v3-sorted-trimmed.bam", checkIfExists: true) ]
    def gff = file("dummy_file.txt")
    def use_gff = false

    QUALIMAP_BAMQC ( input, gff, use_gff )
}
