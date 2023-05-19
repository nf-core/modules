#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPSIFT_DBNSFP } from '../../../../../modules/nf-core/snpsift/dbnsfp/main.nf'

workflow test_snpsift_dbnsfp {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome21_indels_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome21_indels_vcf_gz_tbi'], checkIfExists: true)
    ]
    
    input_dbs = [
        [ id:'databases'],
        file(params.test_data['homo_sapiens']['genome']['dbNSFP_4_1a_21_hg38_txt_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['dbNSFP_4_1a_21_hg38_txt_tbi'], checkIfExists: true)
    ]

    SNPSIFT_DBNSFP ( input, input_dbs )
}

workflow test_snpsift_dbnsfp_uncompressed {
    input = [
        [ id:'tester', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome21_indels_vcf_gz'], checkIfExists: true),
        []
    ]
    
    input_dbs = [
        [id:'databases'],
        file(params.test_data['homo_sapiens']['genome']['dbNSFP_4_1a_21_hg38_txt_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['dbNSFP_4_1a_21_hg38_txt_tbi'], checkIfExists: true)
    ]
    SNPSIFT_DBNSFP ( input, input_dbs )
}
