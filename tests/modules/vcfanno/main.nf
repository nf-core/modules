#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFANNO } from '../../../modules/vcfanno/main.nf'

workflow test_vcfanno {
    
    input = [ 
        [ id:'test_compressed', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
    ]
    
    input_2 = [ [ id:'test_compressed', single_end:false ], // meta map 
                [] ]

    toml = file(params.test_data['homo_sapiens']['genome']['vcfanno_toml'], checkIfExists: true)
    resource_dir = file(params.test_data['homo_sapiens']['genome']['vcfanno_resource_dir'], type:'dir', checkIfExists: true)

    VCFANNO ( input, input_2, toml, resource_dir )
}

workflow test_vcfanno_uncompressed {

    input = [ [ id:'test_uncompressed', single_end:false ], // meta map
                [] ,[] ]
    
    input_2 = [ 
        [ id:'test_uncompressed', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    toml = file(params.test_data['homo_sapiens']['genome']['vcfanno_toml'], checkIfExists: true)
    resource_dir = file(params.test_data['homo_sapiens']['genome']['vcfanno_resource_dir'], type:'dir', checkIfExists: true)

    VCFANNO ( input, input_2, toml, resource_dir )
}