#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFANNO } from '../../../../modules/nf-core/vcfanno/main.nf'

workflow test_vcfanno {

    input = [
        [ id:'test_compressed', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists:true)
    ]
    lua = []
    toml = file(params.test_data['homo_sapiens']['genome']['vcfanno_toml'], checkIfExists: true)
    resources = [
        file("https://github.com/brentp/vcfanno/raw/master/example/exac.vcf.gz", checkIfExists: true),
        file("https://github.com/brentp/vcfanno/raw/master/example/exac.vcf.gz.tbi",checkIfExists: true)
    ]
    VCFANNO ( input, toml, lua, resources )
}

workflow test_vcfanno_uncompressed {

    input = [
        [ id:'test_uncompressed', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        [],
        file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists:true)
    ]
    lua = []
    toml = file(params.test_data['homo_sapiens']['genome']['vcfanno_toml'], checkIfExists: true)
    resources = [
        file("https://github.com/brentp/vcfanno/raw/master/example/exac.vcf.gz", checkIfExists: true),
        file("https://github.com/brentp/vcfanno/raw/master/example/exac.vcf.gz.tbi",checkIfExists: true)
    ]
    VCFANNO ( input, toml, lua, resources )
}
