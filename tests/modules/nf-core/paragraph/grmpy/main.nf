#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARAGRAPH_GRMPY } from '../../../../../modules/nf-core/paragraph/grmpy/main.nf'
include { DELLY_CALL      } from '../../../../../modules/nf-core/delly/call/main.nf'

workflow test_paragraph_grmpy {

    input = Channel.of([
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_genome21_indels_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome21_indels_vcf_gz_tbi'], checkIfExists: true)
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
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    fasta_fai = [
        [ id:'fasta_fai' ],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]

    PARAGRAPH_GRMPY (
        input.combine(bam).combine(manifest),
        fasta,
        fasta_fai
    )
}
