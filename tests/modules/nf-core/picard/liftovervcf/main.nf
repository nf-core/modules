#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_LIFTOVERVCF } from '../../../../../modules/nf-core/picard/liftovervcf/main.nf'

workflow test_picard_liftovervcf {

    input_vcf = [ [ id:'test' ],
                file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
                ]
    dict      = [ [ id:'genome' ],
                file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
                ]
    chain     = [ [ id:'genome' ],
                file(params.test_data['homo_sapiens']['genome']['genome_chain_gz'], checkIfExists: true)
                ]
    fasta     = [ [ id:'genome' ],
                file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
                ]
    PICARD_LIFTOVERVCF ( input_vcf, dict, fasta, chain )
}

workflow test_picard_liftovervcf_stubs {

    input_vcf = [ [ id:'test' ],
                file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
                ]
    dict      = [ [ id:'genome' ],
                file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
                ]
    chain     = [ [ id:'genome' ],
                file(params.test_data['homo_sapiens']['genome']['genome_chain_gz'], checkIfExists: true)
                ]
    fasta     = [ [ id:'genome' ],
                file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
                ]

    PICARD_LIFTOVERVCF ( input_vcf, dict, fasta, chain )
}
