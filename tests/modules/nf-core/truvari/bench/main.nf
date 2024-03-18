#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_BGZIP    } from '../../../../../modules/nf-core/tabix/bgzip/main.nf'
include { TRUVARI_BENCH  } from '../../../../../modules/nf-core/truvari/bench/main.nf'

workflow test_truvari_bench {

    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['simulated_sv'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['simulated_sv_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['simulated_sv2'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['simulated_sv2_tbi'], checkIfExists: true),
        []
    ])

    fasta = [
        [id:'fasta'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [
        [id:'fai'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    TRUVARI_BENCH (
        input,
        fasta,
        fai
    )
}

workflow test_truvari_bench_bed {

    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['simulated_sv'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['simulated_sv_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['simulated_sv2'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['simulated_sv2_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ])

    fasta = [
        [id:'fasta'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]
    fai = [
        [id:'fai'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]

    TRUVARI_BENCH (
        input,
        fasta,
        fai
    )
}
