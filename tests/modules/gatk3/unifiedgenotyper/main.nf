#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK3_UNIFIEDGENOTYPER } from '../../../../modules/gatk3/unifiedgenotyper/main.nf' addParams( options: [:] )

workflow test_gatk3_unifiedgenotyper {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]
    ref = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK3_UNIFIEDGENOTYPER ( input, ref, fai, dict )
}
