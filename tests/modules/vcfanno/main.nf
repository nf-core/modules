#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../modules/untar/main.nf'
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
    resource_dir = [[], file(params.test_data['homo_sapiens']['genome']['vcfanno_tar_gz'], checkIfExists: true) ]

    UNTAR ( resource_dir )
    VCFANNO ( input, input_2, toml, UNTAR.out.untar.map{ it[1] } )
}

workflow test_vcfanno_uncompressed {

    input = [ [ id:'test_uncompressed', single_end:false ], // meta map
                [] ,[] ]

    input_2 = [
        [ id:'test_uncompressed', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    toml = file(params.test_data['homo_sapiens']['genome']['vcfanno_toml'], checkIfExists: true)
    resource_dir = [[], file(params.test_data['homo_sapiens']['genome']['vcfanno_tar_gz'], checkIfExists: true) ]

    UNTAR ( resource_dir )
    VCFANNO ( input, input_2, toml, UNTAR.out.untar.map{ it[1] } )
}
