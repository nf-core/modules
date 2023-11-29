#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CNNSCOREVARIANTS } from '../../../../../modules/nf-core/gatk4/cnnscorevariants/main.nf'

workflow test_gatk4_cnnscorevariants {

    input     = [ [ id:'test' ], // meta map
                    file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true),
                    [],
                    []
                ]
    fasta  = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai    = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict   = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    GATK4_CNNSCOREVARIANTS ( input, fasta, fai, dict, [], [] )
}
