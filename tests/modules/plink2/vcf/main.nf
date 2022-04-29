#!/usr/bin/env nextflow



include { PLINK2_VCF } from '../../../../modules/plink2/vcf/main.nf'

workflow test_plink2_vcf {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ]

    PLINK2_VCF ( input )
}
