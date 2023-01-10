#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAPPY_HAPPY } from '../../../../../modules/nf-core/happy/happy/main.nf'

workflow test_happy_vcf {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_vcf'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true),
        []
    ]

    fasta = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fasta_fai = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    HAPPY_HAPPY (
        input,
        fasta,
        fasta_fai,
        [[],[]],
        [[],[]],
        [[],[]]
    )
}

workflow test_happy_gvcf {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_vcf'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
        [],
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    fasta = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fasta_fai = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    HAPPY_HAPPY (
        input,
        fasta,
        fasta_fai,
        [[],[]],
        [[],[]],
        [[],[]]
    )
}

workflow test_happy_false_positives {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_vcf'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        [],
        []
    ]

    fasta = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fasta_fai = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    false_positives = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    HAPPY_HAPPY (
        input,
        fasta,
        fasta_fai,
        false_positives,
        [[],[]],
        [[],[]]
    )
}

workflow test_happy_stratification {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_vcf'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        [],
        []
    ]

    fasta = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fasta_fai = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    stratification = Channel.of("exon\tgenome.bed").collectFile(name:"stratification.tsv")
        .map{[[id:'test'], it]}
    stratification_beds = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    HAPPY_HAPPY (
        input,
        fasta,
        fasta_fai,
        [[],[]],
        stratification,
        stratification_beds
    )

}
