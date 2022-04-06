#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_ESTIMATELIBRARYCOMPLEXITY } from '../../../../modules/gatk4/estimatelibrarycomplexity/main.nf'

workflow test_gatk4_estimatelibrarycomplexity {

    input = [   [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true)
            ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict  = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_ESTIMATELIBRARYCOMPLEXITY ( input, fasta, fai, dict )
}
