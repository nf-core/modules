#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_PLUGINSPLIT } from '../../../../../modules/nf-core/bcftools/pluginsplit/main.nf'

workflow test_bcftools_pluginsplit {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz_tbi'], checkIfExists: true)
    ]

    samples = Channel.of("normal\t-\tnormal", "tumour\t-\ttumour")
        .collectFile(name:"samples.txt", newLine:true)

    BCFTOOLS_PLUGINSPLIT (
        input,
        samples,
        [],
        [],
        []
    )
}

workflow test_bcftools_pluginsplit_full {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz_tbi'], checkIfExists: true)
    ]

    groups = Channel.of("normal\t-\tnormal", "tumour\t-\ttumour")
        .collectFile(name:"samples.txt", newLine:true)

    regions = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    targets = file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)

    BCFTOOLS_PLUGINSPLIT (
        input,
        [],
        groups,
        regions,
        targets
    )
}
