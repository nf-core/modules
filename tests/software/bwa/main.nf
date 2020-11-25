#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../software/bwa/index/main.nf' addParams( options: [:] )
include { BWA_MEM as BWA_MEM_SE } from '../../../software/bwa/mem/main.nf'  addParams( options: [ publish_dir:'test_single_end' ] )
include { BWA_MEM as BWA_MEM_PE } from '../../../software/bwa/mem/main.nf'  addParams( options: [ publish_dir:'test_paired_end' ] )

workflow test_bwa_index {
    BWA_INDEX ( file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true) )
}

/*
 * Test with single-end data
 */
workflow test_bwa_mem_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/dna/Ecoli_DNA_R1.fastq.gz", checkIfExists: true) ] ]

    BWA_MEM_SE (
        input,
        file("${launchDir}/tests/data/index/E_coli/bwa/NC_010473.fa.{amb,ann,bwt,pac,sa}", checkIfExists: true),
        file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    )
}

/*
 * Test with paired-end data
 */
workflow test_bwa_mem_paired_end {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/fastq/dna/Ecoli_DNA_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/fastq/dna/Ecoli_DNA_R2.fastq.gz", checkIfExists: true) ] ]

    BWA_MEM_PE (
        input,
        file("${launchDir}/tests/data/index/E_coli/bwa/NC_010473.fa.{amb,ann,bwt,pac,sa}", checkIfExists: true),
        file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    )
}
