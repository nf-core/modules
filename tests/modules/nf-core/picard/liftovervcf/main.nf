#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_LIFTOVERVCF } from '../../../../modules/picard/liftovervcf/main.nf'

workflow test_picard_liftovervcf {

    input_vcf = [ [ id:'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
            ]
    dict =  file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    chain =  file(params.test_data['homo_sapiens']['genome']['genome_chain_gz'], checkIfExists: true)
    fasta  = [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]

    PICARD_LIFTOVERVCF ( input_vcf, dict, chain, fasta )
}

workflow test_picard_liftovervcf_stubs {

    input_vcf = [ [ id:'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
            ]
    dict =  file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    chain =  file(params.test_data['homo_sapiens']['genome']['genome_chain_gz'], checkIfExists: true)
    fasta  = [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]

    PICARD_LIFTOVERVCF ( input_vcf, dict, chain, fasta )
}
