#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
def options_align = [args: '--readFilesCommand zcat']
def options_gg = [args: '--genomeSAindexNbases 9']
include { STAR_ALIGN } from '../../../software/star/align/main.nf' addParams( options: options_align )
include { STAR_GENOMEGENERATE } from '../../../software/star/genomegenerate/main.nf'  addParams( options: options_gg )

workflow test_star_genomegenerate {
    fasta = file("${launchDir}/tests/data/fasta/E_coli/GCF_000019425.1_ASM1942v1_genomic.fna", checkIfExists: true)
    gtf = file("${launchDir}/tests/data/gff/GCF_000019425.1_ASM1942v1_genomic.gtf", checkIfExists: true)
    STAR_GENOMEGENERATE ( fasta, gtf )
}

workflow test_star_alignment_single_end {
    fasta = file("${launchDir}/tests/data/fasta/E_coli/GCF_000019425.1_ASM1942v1_genomic.fna", checkIfExists: true)
    gtf = file("${launchDir}/tests/data/gff/GCF_000019425.1_ASM1942v1_genomic.gtf", checkIfExists: true)
    STAR_GENOMEGENERATE ( fasta, gtf )

    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/rna/test_single_end.fastq.gz", checkIfExists: true) ] ]

    STAR_ALIGN( input, STAR_GENOMEGENERATE.out.index, gtf)
}

workflow test_star_alignment_paired_end {
    fasta = file("${launchDir}/tests/data/fasta/E_coli/GCF_000019425.1_ASM1942v1_genomic.fna", checkIfExists: true)
    gtf = file("${launchDir}/tests/data/gff/GCF_000019425.1_ASM1942v1_genomic.gtf", checkIfExists: true)
    STAR_GENOMEGENERATE ( fasta, gtf )

    input = [ [ id:'test', single_end:false ], // meta map
               [ file("${launchDir}/tests/data/fastq/rna/test_R1.fastq.gz", checkIfExists: true),
                 file("${launchDir}/tests/data/fastq/rna/test_R2.fastq.gz", checkIfExists: true) ] ]

    STAR_ALIGN( input, STAR_GENOMEGENERATE.out.index, gtf)
}
