#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLINK_VCF } from '../../../../modules/plink/vcf/main.nf'

workflow test_plink_vcf {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]

    PLINK_VCF ( input )
}
