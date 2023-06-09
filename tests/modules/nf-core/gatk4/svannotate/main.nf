#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SVANNOTATE } from '../../../../../modules/nf-core/gatk4/svannotate/main.nf'
include { MANTA_GERMLINE   } from '../../../../../modules/nf-core/manta/germline/main.nf'

workflow test_gatk4_svannotate {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
        [],
        []
    ]

    fasta     = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)])
    fasta_fai = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)])
    dict      = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)])

    MANTA_GERMLINE(
        input,
        fasta,
        fasta_fai
    )

    GATK4_SVANNOTATE (
        MANTA_GERMLINE.out.diploid_sv_vcf.combine(MANTA_GERMLINE.out.diploid_sv_vcf_tbi, by:0).map({ meta, vcf, tbi -> [ meta, vcf, tbi, [] ]}),
        [],
        [],
        []
    )
}

workflow test_gatk4_svannotate_fasta {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
        [],
        []
    ]

    fasta     = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)])
    fasta_fai = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)])
    dict      = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)])

    MANTA_GERMLINE(
        input,
        fasta,
        fasta_fai
    )

    GATK4_SVANNOTATE (
        MANTA_GERMLINE.out.diploid_sv_vcf.combine(MANTA_GERMLINE.out.diploid_sv_vcf_tbi, by:0).map({ meta, vcf, tbi -> [ meta, vcf, tbi, [] ]}),
        fasta,
        fasta_fai,
        dict
    )
}

workflow test_gatk4_svannotate_bed {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
        [],
        []
    ]

    bed = Channel.of([
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ])

    fasta     = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)])
    fasta_fai = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)])
    dict      = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)])

    MANTA_GERMLINE(
        input,
        fasta,
        fasta_fai
    )

    GATK4_SVANNOTATE (
        MANTA_GERMLINE.out.diploid_sv_vcf
                        .combine(MANTA_GERMLINE.out.diploid_sv_vcf_tbi, by:0)
                        .combine(bed, by:0),
        [],
        [],
        []
    )
}
