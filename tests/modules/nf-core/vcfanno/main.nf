#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from "$moduleDir/modules/nf-core/untar/main.nf"
include { VCFANNO } from "$moduleDir/modules/nf-core/vcfanno/main.nf"

workflow test_vcfanno {

    input = [
        [ id:'test_compressed', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
    ]

    toml = file(params.test_data['homo_sapiens']['genome']['vcfanno_toml'], checkIfExists: true)
    resource_dir = [[], file(params.test_data['homo_sapiens']['genome']['vcfanno_tar_gz'], checkIfExists: true) ]

    UNTAR ( resource_dir )
    VCFANNO ( input, toml, UNTAR.out.untar.map{ it[1] } )
}

workflow test_vcfanno_uncompressed {

    input = [
        [ id:'test_uncompressed', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ]

    toml = file(params.test_data['homo_sapiens']['genome']['vcfanno_toml'], checkIfExists: true)
    resource_dir = [[], file(params.test_data['homo_sapiens']['genome']['vcfanno_tar_gz'], checkIfExists: true) ]

    UNTAR ( resource_dir )
    VCFANNO ( input, toml, UNTAR.out.untar.map{ it[1] } )
}
