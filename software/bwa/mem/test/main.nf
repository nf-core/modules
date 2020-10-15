#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_MEM as BWA_MEM_SE } from '../main.nf'  addParams( options: [ publish_dir:'test_single_end' ] )
include { BWA_MEM as BWA_MEM_PE } from '../main.nf'  addParams( options: [ publish_dir:'test_paired_end' ] )

/*
 * Test with single-end data
 */
workflow test_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${baseDir}/input/Ecoli_DNA_R1.fastq.gz", checkIfExists: true) ] ]

    BWA_MEM_SE (
        input,
        file("${baseDir}/input/index/NC_010473.fa.{amb,ann,bwt,pac,sa}", checkIfExists: true),
        file("${baseDir}/input/NC_010473.fa", checkIfExists: true)
    )
}

/*
 * Test with paired-end data
 */
workflow test_paired_end {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${baseDir}/input/Ecoli_DNA_R1.fastq.gz", checkIfExists: true),
                file("${baseDir}/input/Ecoli_DNA_R2.fastq.gz", checkIfExists: true) ] ]

    BWA_MEM_PE (
        input,
        file("${baseDir}/input/index/NC_010473.fa.{amb,ann,bwt,pac,sa}", checkIfExists: true),
        file("${baseDir}/input/NC_010473.fa", checkIfExists: true)
    )
}

workflow {
    test_single_end()
    test_paired_end()
}
