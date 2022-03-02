#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_BAM2FQ } from '../../../../modules/samtools/bam2fq/main.nf'

workflow test_samtools_bam2fq_nosplit {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_umi_converted_bam'], checkIfExists: true) ]
    split = false

    SAMTOOLS_BAM2FQ ( input, split )
}


workflow test_samtools_bam2fq_withsplit {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_umi_converted_bam'], checkIfExists: true) ]
    split = true

    SAMTOOLS_BAM2FQ ( input, split )
}
