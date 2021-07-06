#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IVAR_VARIANTS } from '../../../../software/ivar/variants/main.nf' addParams([:])

workflow test_ivar_variants_no_gff_no_mpileup {
    params.gff          = false
    params.save_mpileup = false

    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) 
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    dummy = file("dummy_file.txt")
    
    IVAR_VARIANTS ( input, fasta, dummy )
}

workflow test_ivar_variants_no_gff_with_mpileup {
    params.gff          = false
    params.save_mpileup = true

    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) 
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    dummy = file("dummy_file.txt")

    IVAR_VARIANTS ( input, fasta, dummy )
}

workflow test_ivar_variants_with_gff_with_mpileup {
    params.gff          = true
    params.save_mpileup = true
    
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) 
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    gff   = file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true)
    
    IVAR_VARIANTS ( input, fasta, gff )
}
