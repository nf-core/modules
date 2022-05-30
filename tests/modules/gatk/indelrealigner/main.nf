#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_INDELREALIGNER } from '../../../../modules/gatk/indelrealigner/main.nf'

// TODO add REalignerTargetCrator


workflow test_gatk_indelrealigner {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bai'], checkIfExists: true),
        GATK_REALIGNERTARGETCREATOR.out.intervals
    ]

    reference = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GATK_INDELREALIGNER ( input, reference, [] )
}
