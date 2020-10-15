#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../main.nf' addParams( options: [:] )

workflow test {
    BWA_INDEX ( file("${baseDir}/input/NC_010473.fa", checkIfExists: true) )
}

workflow {
    test()
}
