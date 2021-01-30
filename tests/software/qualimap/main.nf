#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUALIMAP_BAMQC } from '../../../software/qualimap/bamqc/main.nf' addParams( options: [:] )

workflow test_qualimap_bamqc {

    def input = []
    def gff = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ]
    gff = [ file("dummy_file.txt") ]

    QUALIMAP_BAMQC ( input, gff )
}
