#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_ALIGN } from '../../../software/star/align/main.nf' addParams( options: [:] )
include { STAR_GENOMEGENERATE } from '../../../software/star/genomegenerate/main.nf'  addParams( options: [:] )

workflow test_star_genomegenerate {
    fasta = file("${launchDir}/tests/data/fasta/E_coli/GCF_000019425.1_ASM1942v1_genomic.fna", checkIfExists: true)
    gtf = file("${launchDir}/tests/data/gff/GCF_000019425.1_ASM1942v1_genomic.gtf", checkIfExists: true)
    STAR_GENOMEGENERATE ( fasta, gtf )
}

// workflow test_star_alignment_single_end {

//     fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
//     BOWTIE_INDEX ( fasta )

//     def input = []
//     input = [ [ id:'test', single_end:true ], // meta map
//               [ file("${launchDir}/tests/data/fastq/rna/test_R1.fastq.gz", checkIfExists: true) ] ]
//     BOWTIE_ALIGN ( input, BOWTIE_INDEX.out.index )
// }

// workflow test_star_alignment_paired_end {
//     fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
//     BOWTIE_INDEX ( fasta )

//     def input = []
//     input = [ [ id:'test', single_end:false ], // meta map
//               [ file("${launchDir}/tests/data/fastq/rna/test_R1.fastq.gz", checkIfExists: true),
//                 file("${launchDir}/tests/data/fastq/rna/test_R2.fastq.gz", checkIfExists: true) ] ]
//     BOWTIE_ALIGN ( input, BOWTIE_INDEX.out.index )
// }
