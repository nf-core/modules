#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RHOCALL_ANNOTATE } from '../../../../../modules/nf-core/rhocall/annotate/main.nf'
include { BCFTOOLS_ROH } from '../../../../../modules/nf-core/bcftools/roh/main.nf'

workflow test_rhocall_annotate {

    input = [ [ id:'test' ], // meta map
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)]

    af_file = [[],[]]
    gen_map = []
    regions = []
    targets = []
    samples = []

    BCFTOOLS_ROH ( input, af_file, gen_map, regions, samples, targets )
    RHOCALL_ANNOTATE ( input, BCFTOOLS_ROH.out.roh, [])

}

workflow test_rhocall_annotate_stub {

    input = [ [ id:'test' ], // meta map
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)]

    af_file = [[],[]]
    gen_map = []
    regions = []
    targets = []
    samples = []

    BCFTOOLS_ROH ( input, af_file, gen_map, regions, samples, targets )
    RHOCALL_ANNOTATE ( input, BCFTOOLS_ROH.out.roh, [])

}
