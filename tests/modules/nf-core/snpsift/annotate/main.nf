#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPSIFT_ANNOTATE } from '../../../../../modules/nf-core/snpsift/annotate/main.nf'

workflow test_snpsift_annotate {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
    ]
    
    input_dbs = [
        [ id:'databases'],
        file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    ]
    
    SNPSIFT_ANNOTATE ( input, input_dbs )
}

workflow test_snpsift_annotate_uncompressed{
    input = [
        [ id:'tester', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ]
    
    input_dbs = [
        [id:'databases'],
        file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    ]
    
    SNPSIFT_ANNOTATE ( input, input_dbs )
}