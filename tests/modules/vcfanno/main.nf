#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFANNO } from '../../../modules/vcfanno/main.nf'

workflow test_vcfanno {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
    ]

    resource_dir = file(params.test_data['homo_sapiens']['genome']['vcfanno_resource_dir'], type:'dir', checkIfExists: true)

    VCFANNO ( input, resource_dir )
}
