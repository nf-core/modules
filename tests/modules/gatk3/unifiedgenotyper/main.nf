#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK3_UNIFIEDGENOTYPER } from '../../../../modules/gatk3/unifiedgenotyper/main.nf' addParams( options: [:] )

workflow test_gatk3_unifiedgenotyper {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
    ref = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GATK3_UNIFIEDGENOTYPER ( input, ref )
}
