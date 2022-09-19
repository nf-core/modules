#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVDB_QUERY } from '../../../../modules/svdb/query/main.nf'

workflow test_svdb_query {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_sv_vcf'], checkIfExists: true) ]
            ]

    vcf_db = [
                file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_sv_vcf_gz'], checkIfExists: true)
            ]

    in_occs = ['AC']
    in_frqs = ['AF']
    out_occs = ['gnomad_svAC']
    out_frqs = ['gnomad_svAF']

    SVDB_QUERY ( input, in_occs, in_frqs, out_occs, out_frqs, vcf_db )
}

workflow test_svdb_query_multiple {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_sv_vcf'], checkIfExists: true) ]
            ]

    vcf_db = [
                file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_sv_vcf_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['genome']['gnomad2_r2_1_1_sv_vcf_gz'], checkIfExists: true)
            ]

    in_occs = ['AC','AC']
    in_frqs = ['AF','AF']
    out_occs = ['gnomad_svAC','gnomad_svAC']
    out_frqs = ['gnomad_svAF','gnomad_svAF']

    SVDB_QUERY ( input, in_occs, in_frqs, out_occs, out_frqs, vcf_db )
}
