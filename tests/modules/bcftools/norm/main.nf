#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_NORM } from '../../../../modules/bcftools/norm/main.nf'

workflow test_bcftools_norm {
    
    input = [ [ id:'test2', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BCFTOOLS_NORM ( input, fasta )
}
