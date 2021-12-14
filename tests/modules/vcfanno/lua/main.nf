#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFANNO_LUA } from '../../../../modules/vcfanno/lua/main.nf'

workflow test_vcfanno_lua {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    vcfanno_config    = file(params.test_data['sarscov2']['illumina']['test_vcfanno_conf_toml'], checkIfExists: true)
    vcfanno_functions = []

    VCFANNO_LUA ( input, vcfanno_config, vcfanno_functions )
}
