#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_VALIDATE_SMALL_VARIANTS } from '../../../../subworkflows/nf-core/vcf_validate_small_variants/main.nf'
include { UNTAR             } from '../../../../modules/nf-core/untar/main.nf'

workflow test_vcf_validate_small_variants_happy {

    input = Channel.of([
        [id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz_tbi'], checkIfExists: true),
    ])

    beds = Channel.of([
        [id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true),
        []
    ])

    fasta = Channel.of([
        [id:'fasta'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ])
    fasta_fai = Channel.of([
        [id:'fasta_fai'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ])

    VCF_VALIDATE_SMALL_VARIANTS (
        input,
        beds,
        fasta,
        fasta_fai,
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]],
        "happy"
    )
}

workflow test_vcf_validate_small_variants_vcfeval {

    input = Channel.of([
        [id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz_tbi'], checkIfExists: true),
    ])

    beds = Channel.of([
        [id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true),
        []
    ])

    compressed_sdf = [
        [ id:'sdf' ],
        file(params.test_data['homo_sapiens']['genome']['genome_21_sdf'])
    ]

    sdf = UNTAR( compressed_sdf ).untar

    VCF_VALIDATE_SMALL_VARIANTS (
        input,
        beds,
        [[],[]],
        [[],[]],
        sdf,
        [[],[]],
        [[],[]],
        [[],[]],
        "vcfeval"
    )
}

workflow test_vcf_validate_small_variants_all {

    input = Channel.of([
        [id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz_tbi'], checkIfExists: true),
    ])

    beds = Channel.of([
        [id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true),
        []
    ])

    fasta = Channel.of([
        [id:'fasta'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ])
    fasta_fai = Channel.of([
        [id:'fasta_fai'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ])

    compressed_sdf = [
        [ id:'sdf' ],
        file(params.test_data['homo_sapiens']['genome']['genome_21_sdf'])
    ]

    sdf = UNTAR( compressed_sdf ).untar

    VCF_VALIDATE_SMALL_VARIANTS (
        input,
        beds,
        fasta,
        fasta_fai,
        sdf,
        [[],[]],
        [[],[]],
        [[],[]],
        "happy,vcfeval"
    )
}

workflow test_vcf_validate_small_variants_happy_optional {

    input = Channel.of([
        [id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz_tbi'], checkIfExists: true),
    ])

    beds = Channel.of([
        [id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true),
        []
    ])

    fasta = Channel.of([
        [id:'fasta'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ])
    fasta_fai = Channel.of([
        [id:'fasta_fai'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ])

    stratification = Channel.of("exon\tgenome.bed").collectFile(name:"stratification.tsv")
        .map{[[id:'test'], it]}
    stratification_beds = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    false_positives = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_antitarget_bed'], checkIfExists: true)
    ]

    VCF_VALIDATE_SMALL_VARIANTS (
        input,
        beds,
        fasta,
        fasta_fai,
        [[],[]],
        false_positives,
        stratification,
        stratification_beds,
        "happy"
    )
}
