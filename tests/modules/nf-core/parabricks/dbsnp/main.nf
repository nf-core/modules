#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_DBSNP } from '../../../../../modules/nf-core/parabricks/dbsnp/main.nf'

workflow test_parabricks_dbsnp {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_vcf'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)
    ]

    PARABRICKS_DBSNP ( input )
}
