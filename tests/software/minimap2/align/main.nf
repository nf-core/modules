#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIMAP2_ALIGN } from '../../../../software/minimap2/align/main.nf' addParams( options: [:] )

workflow test_minimap2_align_single_end {

    def fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/rna/test_R1.fastq.gz", checkIfExists: true) ] ]
    MINIMAP2_ALIGN ( input, fasta )
}

workflow test_minimap2_align_paired_end {
    
    def fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/fastq/rna/test_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/fastq/rna/test_R2.fastq.gz", checkIfExists: true) ] ]
    MINIMAP2_ALIGN ( input, fasta )
}

workflow test_minimap2_align_pairwise {
    
    def fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ fasta ] ]

    MINIMAP2_ALIGN ( input, fasta )
}
