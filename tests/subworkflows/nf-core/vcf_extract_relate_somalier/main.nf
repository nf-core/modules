#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_EXTRACT_RELATE_SOMALIER } from '../../../../subworkflows/nf-core/vcf_extract_relate_somalier/main.nf'

workflow test_vcf_extract_relate_somalier_minimum {

    vcfs = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        [],
        1
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    somalier_sites = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/sites_chr21.hg38.vcf.gz", checkIfExists: true)
    peds = Channel.of([
        [id:"test"],
        []
    ])
    sample_groups = []

    VCF_EXTRACT_RELATE_SOMALIER (
        vcfs,
        fasta,
        fasta_fai,
        somalier_sites,
        peds,
        sample_groups
    )
}

workflow test_vcf_extract_relate_somalier_index {

    vcfs = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
        1
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    somalier_sites = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/sites_chr21.hg38.vcf.gz", checkIfExists: true)
    peds = Channel.of([
        [id:"test"],
        []
    ])

    sample_groups = []

    VCF_EXTRACT_RELATE_SOMALIER (
        vcfs,
        fasta,
        fasta_fai,
        somalier_sites,
        peds,
        sample_groups
    )
}

workflow test_vcf_extract_relate_somalier_ped {

    vcfs = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
        1
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    somalier_sites = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/sites_chr21.hg38.vcf.gz", checkIfExists: true)
    peds = Channel.of([
        [id:"test"],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/family.ped", checkIfExists: true)
    ])
    sample_groups = []

    VCF_EXTRACT_RELATE_SOMALIER (
        vcfs,
        fasta,
        fasta_fai,
        somalier_sites,
        peds,
        sample_groups
    )
}

workflow test_vcf_extract_relate_somalier_mixed_combine {

    vcfs = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        [],
        2
    ],[
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true),
        2
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    somalier_sites = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/sites_chr21.hg38.vcf.gz", checkIfExists: true)
    peds = Channel.of([
        [id:"test"],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/family.ped", checkIfExists: true)
    ])
    sample_groups = Channel.of("disease_103,testN").collectFile(name:"sample_groups.txt")

    VCF_EXTRACT_RELATE_SOMALIER (
        vcfs,
        fasta,
        fasta_fai,
        somalier_sites,
        peds,
        sample_groups
    )
}

workflow test_vcf_extract_relate_somalier_mixed_no_combine {

    vcfs = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
        1
    ],[
        [id:"test2"],
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        [],
        1
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    somalier_sites = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/sites_chr21.hg38.vcf.gz", checkIfExists: true)
    peds = Channel.of([
        [id:"test"],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/family.ped", checkIfExists: true)
    ],
    [
        [id:"test2"],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/family.ped", checkIfExists: true)
    ])
    sample_groups = []

    VCF_EXTRACT_RELATE_SOMALIER (
        vcfs,
        fasta,
        fasta_fai,
        somalier_sites,
        peds,
        sample_groups
    )
}

workflow test_vcf_extract_relate_somalier_joint_vcf {

    vcfs = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz_tbi'], checkIfExists: true),
        []
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    somalier_sites = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/sites_chr21.hg38.vcf.gz", checkIfExists: true)
    peds = Channel.of([
        [id:"test"],
        []
    ])
    sample_groups = []

    VCF_EXTRACT_RELATE_SOMALIER (
        vcfs,
        fasta,
        fasta_fai,
        somalier_sites,
        peds,
        sample_groups
    )
}

workflow test_vcf_extract_relate_somalier_mixed_combine_no_count {

    vcfs = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        [],
        []
    ],[
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true),
        []
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    somalier_sites = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/sites_chr21.hg38.vcf.gz", checkIfExists: true)
    peds = Channel.of([
        [id:"test"],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/family.ped", checkIfExists: true)
    ])
    sample_groups = Channel.of("disease_103,testN").collectFile(name:"sample_groups.txt")

    VCF_EXTRACT_RELATE_SOMALIER (
        vcfs,
        fasta,
        fasta_fai,
        somalier_sites,
        peds,
        sample_groups
    )
}
