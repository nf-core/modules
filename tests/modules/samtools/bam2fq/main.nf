#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_BAM2FQ } from '../../../../modules/samtools/bam2fq/main.nf' addParams( options: [:] )

workflow test_samtools_bam2fq {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    SAMTOOLS_BAM2FQ ( input )
}
