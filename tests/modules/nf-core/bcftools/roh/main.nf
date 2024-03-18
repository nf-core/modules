#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_ROH } from '../../../../../modules/nf-core/bcftools/roh/main.nf'

workflow test_bcftools_roh {

    input = [ [ id:'test' ], // meta map
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)]

    af_file = [[],[]]
    gen_map = []
    regions = []
    targets = []
    samples = []

    BCFTOOLS_ROH ( input, af_file, gen_map, regions, samples, targets )
}

workflow test_bcftools_roh_stub {

    input = [ [ id:'test' ], // meta map
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)]

    af_file = [[],[]]
    gen_map = []
    regions = []
    targets = []
    samples = []

    BCFTOOLS_ROH ( input, af_file, gen_map, regions, samples, targets )
}
