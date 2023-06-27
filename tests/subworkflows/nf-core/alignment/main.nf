#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALIGNMENT } from '../../../../subworkflows/nf-core/alignment/main.nf'

workflow test_alignment {
    
    // channels enable parralle: https://www.nextflow.io/docs/latest/faq.html?highlight=parallel
    
    fastqs = [
    [[id:'gene', single_end:false], [params.test_data['uncollapsed_bam_generation']['merged_fastq']['merged_1'], params.test_data['uncollapsed_bam_generation']['merged_fastq']['merged_2']]]
    ]
    reference = [
        [id:'reference'], 
        file('test_nucleo/reference/')
    ]
    fastqs = ch_fastq = Channel.fromList(fastqs)
    ALIGNMENT ( fastqs, reference, 1)
}
