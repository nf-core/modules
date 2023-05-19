#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//keep --no-verson argument, otherwise md5 will change on each execution
include { BCFTOOLS_MERGE                        } from '../../../../../modules/nf-core/bcftools/merge/main.nf'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_GVCF } from '../../../../../modules/nf-core/bcftools/merge/main.nf'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_BCF  } from '../../../../../modules/nf-core/bcftools/merge/main.nf'

workflow test_bcftools_merge {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz'], checkIfExists: true) ],
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz_tbi'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz_tbi'], checkIfExists: true) ]
            ]

    bed = []
    fasta = [[],[]]
    fai = [[],[]]

    BCFTOOLS_MERGE ( input, fasta, fai, bed )
}

workflow test_bcftools_merge_bed {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz'], checkIfExists: true) ],
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz_tbi'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz_tbi'], checkIfExists: true) ]
            ]

    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    fasta = [[],[]]
    fai = [[],[]]

    BCFTOOLS_MERGE ( input, fasta, fai, bed )
}

workflow test_bcftools_merge_gvcf {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true) ],
              [ file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf_gz_tbi'], checkIfExists: true) ]
            ]

    bed = []
    fasta = [ [ id:'genome' ], // meta map
                 file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]
    fai   = [ [ id:'genome' ], // meta map
                 file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
            ]
    BCFTOOLS_MERGE_GVCF ( input, fasta, fai, bed )
}

workflow test_bcftools_merge_bcf {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz'], checkIfExists: true) ],
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz_tbi'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz_tbi'], checkIfExists: true) ]
            ]

    bed = []
    fasta = [[],[]]
    fai = [[],[]]

    BCFTOOLS_MERGE_BCF ( input, fasta, fai, bed )
}
