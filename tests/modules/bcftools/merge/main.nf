#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//keep --no-verson argument, otherwise md5 will change on each execution
include { BCFTOOLS_MERGE                        } from '../../../../modules/bcftools/merge/main.nf'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_GVCF } from '../../../../modules/bcftools/merge/main.nf'

workflow test_bcftools_merge {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz'], checkIfExists: true) ],
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz_tbi'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz_tbi'], checkIfExists: true) ],
              []
            ]

    fasta = []
    fasta_fai = []

    BCFTOOLS_MERGE ( input, fasta, fasta_fai )
}

workflow test_bcftools_merge_bed {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz'], checkIfExists: true) ],
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz_tbi'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz_tbi'], checkIfExists: true) ],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]

    fasta = []
    fasta_fai = []

    BCFTOOLS_MERGE ( input, fasta, fasta_fai )
}

workflow test_bcftools_merge_gvcf {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true) ],
              [ file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf_gz_tbi'], checkIfExists: true) ],
              []
            ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    BCFTOOLS_MERGE_GVCF ( input, fasta, fasta_fai )
}