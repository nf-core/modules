#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE_INDEX } from '../../../software/bowtie/index/main.nf' addParams( options: [:] )
include { BOWTIE_ALIGN } from '../../../software/bowtie/align/main.nf'  addParams( options: [:] )



workflow test_bowtie_index {
    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    BOWTIE_INDEX ( fasta ) 
}

workflow test_bowtie_alignment_single_end {

    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    BOWTIE_INDEX ( fasta ) 

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/dna/Ecoli_DNA_R1.fastq.gz", checkIfExists: true) ] ]
    BOWTIE_ALIGN ( input, BOWTIE_INDEX.out.index )
}

