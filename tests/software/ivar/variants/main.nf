#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IVAR_VARIANTS } from '../../../../software/ivar/variants/main.nf' addParams([:])

workflow test_ivar_variants_no_gff_no_mpileup {
    params.gff          = false
    params.save_mpileup = false

    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    dummy = file("dummy_file.txt")
    
    IVAR_VARIANTS ( input, fasta, dummy )
}

workflow test_ivar_variants_no_gff_with_mpileup {
    params.gff          = false
    params.save_mpileup = true

    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    dummy = file("dummy_file.txt")

    IVAR_VARIANTS ( input, fasta, dummy )
}

workflow test_ivar_variants_with_gff_with_mpileup {
    params.gff          = true
    params.save_mpileup = true
    
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    gff   = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.gff3", checkIfExists: true)
    
    IVAR_VARIANTS ( input, fasta, gff )
}
