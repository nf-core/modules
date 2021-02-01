#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE_BUILD } from '../../../software/bowtie/build/main.nf' addParams( options: [:] )
include { BOWTIE_ALIGN } from '../../../software/bowtie/align/main.nf'  addParams( options: [:] )

workflow test_bowtie_build {
    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    BOWTIE_BUILD ( fasta )
}

workflow test_bowtie_alignment_single_end {

    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    BOWTIE_BUILD ( fasta )

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/rna/test_R1.fastq.gz", checkIfExists: true) ] ]
    BOWTIE_ALIGN ( input, BOWTIE_BUILD.out.index )
}

workflow test_bowtie_alignment_paired_end {
    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    BOWTIE_BUILD ( fasta )

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/fastq/rna/test_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/fastq/rna/test_R2.fastq.gz", checkIfExists: true) ] ]
    BOWTIE_ALIGN ( input, BOWTIE_BUILD.out.index )
}
