#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTMULTIPLEMETRICS } from '../main.nf' addParams( options: [:] )

workflow test {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${baseDir}/input/test.paired_end.sorted.bam", checkIfExists: true) ]

    PICARD_COLLECTMULTIPLEMETRICS (
        input,
        file("${baseDir}/input/NC_010473.fa", checkIfExists: true)
    )
}

workflow {
    test()
}
