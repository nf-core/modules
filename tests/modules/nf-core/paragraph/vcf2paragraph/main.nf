#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARAGRAPH_VCF2PARAGRAPH } from '../../../../../modules/nf-core/paragraph/vcf2paragraph/main.nf'

workflow test_paragraph_vcf2paragraph {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true)
    ]

    fasta = [[], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]

    PARAGRAPH_VCF2PARAGRAPH (
        input,
        fasta
    )
}
