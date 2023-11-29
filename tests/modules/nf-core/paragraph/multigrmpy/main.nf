#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARAGRAPH_MULTIGRMPY } from '../../../../../modules/nf-core/paragraph/multigrmpy/main.nf'

workflow test_paragraph_multigrmpy {

    input = Channel.of([
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true)
    ])

    manifest = Channel.of(
        "id\tpath\tdepth\tread length\ntest\ttest.paired_end.sorted.cram\t0.73\t150"
    ).collectFile(name:"manifest.txt")

    bam = Channel.of([
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ])

    fasta = [
        [ id:'fasta' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fasta_fai = [
        [ id:'fasta_fai' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    PARAGRAPH_MULTIGRMPY (
        input.combine(bam).combine(manifest),
        fasta,
        fasta_fai
    )
}

workflow test_paragraph_multigrmpy_empty_vcf {

    input = Channel.of([
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['illumina']['empty_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['empty_vcf_gz_tbi'], checkIfExists: true)
    ])

    manifest = Channel.of(
        "id\tpath\tdepth\tread length\ntest\ttest.paired_end.sorted.cram\t0.73\t150"
    ).collectFile(name:"manifest.txt")

    bam = Channel.of([
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ])

    fasta = [
        [ id:'fasta' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fasta_fai = [
        [ id:'fasta_fai' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    PARAGRAPH_MULTIGRMPY (
        input.combine(bam).combine(manifest),
        fasta,
        fasta_fai
    )
}
