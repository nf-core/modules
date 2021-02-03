#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_MEM } from '../../../software/bwa/mem/main.nf' addParams( options: [:] )

/*
 * Test with single-end data
 */
workflow test_bwa_mem_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/dna/Ecoli_DNA_R1.fastq.gz", checkIfExists: true) ] ]

    BWA_MEM (
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

    BWA_MEM (
        input,
        file("${launchDir}/tests/data/index/E_coli/bwa/NC_010473.fa.{amb,ann,bwt,pac,sa}", checkIfExists: true),
        file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    )
}
