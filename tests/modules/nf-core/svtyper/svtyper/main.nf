#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVTYPER_SVTYPER } from '../../../../../modules/nf-core/svtyper/svtyper/main.nf'


workflow test_svtyper_svtyper {
    fa = [
        [id:"ref"],
        fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]
    fai = [
        [id:"ref"],
        fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
    ]

    SVTYPER_SVTYPER ( input, fa, fai )
}
