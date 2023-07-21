#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_UNIFIEDGENOTYPER } from '../../../../modules/nf-core/gatk/unifiedgenotyper/main.nf'
include { TOPAS_GENCONS } from '../../../../../modules/nf-core/topas/gencons/main.nf'

workflow test_topas_gencons {

    input_gatk     = [ [ id:'test' ], // meta map
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

    GATK_UNIFIEDGENOTYPER ( input_gatk, fasta, fai, dict, [[],[]], [[],[]], [[],[]], [[],[]])

    gencons_vcf = GATK_UNIFIEDGENOTYPER.out.vcf
    gencons_fasta =[ [ id:'test' ], // meta map
                     file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
                    ]


    TOPAS_GENCONS ( gencons_vcf, gencons_fasta )
}
