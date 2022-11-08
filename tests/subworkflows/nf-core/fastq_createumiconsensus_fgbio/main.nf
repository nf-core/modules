#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_CREATEUMICONSENSUS_FGBIO } from '../../../../subworkflows/nf-core/fastq_createumiconsensus_fgbio/main.nf'

workflow test_fastq_createumiconsensus_fgbio {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    FASTQ_CREATEUMICONSENSUS_FGBIO ( input )
}
