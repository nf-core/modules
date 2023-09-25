#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_SPLIT } from '../../../../../modules/nf-core/bcftools/split/main.nf'

workflow test_bcftools_split {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true)
    ]

    BCFTOOLS_SPLIT ( input )
}
