#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF2CYTOSURE } from '../../../../modules/nf-core/vcf2cytosure/main.nf'

workflow test_vcf2cytosure {

    sv_vcf = [
        [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['na24385_chr22_sv_vcf'], checkIfExists: true) ]
    ]
    coverage = [
        [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['na24385_chr22_coverage'], checkIfExists: true) ]
    ]
    VCF2CYTOSURE ( sv_vcf, coverage, [[],[]], [[],[]], [] )

}
