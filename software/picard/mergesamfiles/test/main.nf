#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_MERGESAMFILES } from '../main.nf' addParams( options: [:] )

workflow test {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${baseDir}/input/test.paired_end.sorted.bam", checkIfExists: true),
                file("${baseDir}/input/test.paired_end.COPY.sorted.bam", checkIfExists: true), ] ]

    PICARD_MERGESAMFILES ( input )
}

workflow {
    test()
}
