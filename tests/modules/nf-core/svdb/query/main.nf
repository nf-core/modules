#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVDB_QUERY as SVDB_QUERY_SINGLE } from '../../../../../modules/nf-core/svdb/query/main.nf'
include { SVDB_QUERY as SVDB_QUERY_MULTIPLE } from '../../../../../modules/nf-core/svdb/query/main.nf'

workflow test_svdb_query {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_sv_vcf'], checkIfExists: true) ]
            ]

    vcf_db = [
                file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_sv_vcf_gz'], checkIfExists: true)
            ]

    SVDB_QUERY_SINGLE ( input, vcf_db )
}

workflow test_svdb_query_multiple {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_sv_vcf'], checkIfExists: true) ]
            ]

    vcf_db = [
                file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_sv_vcf_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['genome']['gnomad2_r2_1_1_sv_vcf_gz'], checkIfExists: true)
            ]

    SVDB_QUERY_MULTIPLE ( input, vcf_db )
}
