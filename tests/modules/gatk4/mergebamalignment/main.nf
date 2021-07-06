#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MERGEBAMALIGNMENT } from '../../../../software/gatk4/mergebamalignment/main.nf' addParams( options: [:] )

workflow test_gatk4_mergebamalignment {
    input    = [ [ id:'test' ], // meta map
                 file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
               ]
    unmapped = file(params.test_data['sarscov2']['illumina']['test_unaligned_bam'], checkIfExists: true)
    fasta    = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    dict     = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_MERGEBAMALIGNMENT ( input, unmapped, fasta, dict )
}
