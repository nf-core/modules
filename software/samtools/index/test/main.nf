#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from '../main.nf' addParams( options: [:] )

workflow test {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${baseDir}/input/test.paired_end.sorted.bam", checkIfExists: true) ]

    SAMTOOLS_INDEX ( input )
}

workflow {
    test()
}
