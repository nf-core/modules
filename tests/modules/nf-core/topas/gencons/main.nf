#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_UNIFIEDGENOTYPER as GATK_UNIFIEDGENOTYPERSNPS } from '../../../../modules/nf-core/gatk/unifiedgenotyper/main.nf'
include { GATK_UNIFIEDGENOTYPER as GATK_UNIFIEDGENOTYPERINDELS } from '../../../../modules/nf-core/gatk/unifiedgenotyper/main.nf'
include { TOPAS_GENCONS         } from '../../../../../modules/nf-core/topas/gencons/main.nf'

workflow test_topas_gencons {

    input_gatk = [ [ id:'test' ], // meta map
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    dict = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
    ]

    GATK_UNIFIEDGENOTYPERSNPS ( input_gatk, fasta, fai, dict, [[],[]], [[],[]], [[],[]], [[],[]])

    gencons_vcf = GATK_UNIFIEDGENOTYPERSNPS.out.vcf
    gencons_vcf_indels = [[],[]]
    gencons_fasta =[ [ id:'test' ], // meta map
                     file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
                    ]
    gencons_vcf_output = false

    TOPAS_GENCONS ( gencons_vcf, gencons_vcf_indels, gencons_fasta, [[],[]], gencons_vcf_output)
}

workflow test_topas_gencons_fai {

    input_gatk = [ [ id:'test' ], // meta map
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    dict = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
    ]

    GATK_UNIFIEDGENOTYPERSNPS ( input_gatk, fasta, fai, dict, [[],[]], [[],[]], [[],[]], [[],[]])

    gencons_vcf = GATK_UNIFIEDGENOTYPERSNPS.out.vcf
    gencons_vcf_indels = [[],[]]
    gencons_fasta =[ [ id:'test' ], // meta map
                     file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
                    ]
    gencons_vcf_output = false

    TOPAS_GENCONS ( gencons_vcf, gencons_vcf_indels, gencons_fasta, fai, gencons_vcf_output)
}

workflow test_topas_gencons_indels {

    input_gatk = [ [ id:'test' ], // meta map
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    dict = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
    ]

    GATK_UNIFIEDGENOTYPERSNPS ( input_gatk, fasta, fai, dict, [[],[]], [[],[]], [[],[]], [[],[]])
    GATK_UNIFIEDGENOTYPERINDELS ( input_gatk, fasta, fai, dict, [[],[]], [[],[]], [[],[]], [[],[]])


    gencons_vcf = GATK_UNIFIEDGENOTYPERSNPS.out.vcf
    gencons_vcf_indels = GATK_UNIFIEDGENOTYPERINDELS.out.vcf
    gencons_fasta =[ [ id:'test' ], // meta map
                     file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
                    ]
    gencons_vcf_output = true

    TOPAS_GENCONS ( gencons_vcf, gencons_vcf_indels, gencons_fasta, [[],[]], gencons_vcf_output)
}
