#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HISAT2_BUILD } from '../../../software/hisat2/build/main.nf' addParams( options: [:] )

workflow test_hisat2_build {
    
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true) ]

    HISAT2_BUILD ( input )
}
