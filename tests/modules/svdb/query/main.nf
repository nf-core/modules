#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVDB_QUERY } from '../../../../modules/svdb/query/main.nf'

workflow test_svdb_query {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_sv_vcf'], checkIfExists: true) ]
            ]

    allele_freq_count = [
        "AC", "AF", "gnomad_svAC", "gnomad_svAF"
    ]

    vcf_db = [
                file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_sv_vcf_gz'], checkIfExists: true)
            ]

    SVDB_QUERY ( input, allele_freq_count, vcf_db )
}
