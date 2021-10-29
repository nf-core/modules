#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLINK_VCF } from '../../../modules/plink/vcf/main.nf' addParams( options: [:] )
include { ADMIXTURE } from '../../../modules/admixture/main.nf' addParams( options: [:] )

workflow test_admixture {
    
    input = [ [ id:'test', single_end:false ], // meta map
             file(params.test_data['homosapiens']['genome']['mills_and_1000g_indels_vcf_gz'], checkIfExists: true) ]
    PLINK_VCF ( input )

    ADMIXTURE ( PLINK_VCF.out.bed, PLINK_VCF.out.fam )
}
