#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VGAN_HAPLOCART } from '../../../../../modules/nf-core/vgan/haplocart/main.nf'
include { SAMTOOLS_BAM2FQ } from '../../../../../modules/nf-core/samtools/bam2fq/main.nf'


workflow test_vgan_haplocart {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['mitochon_standin_recalibrated_sorted_bam'], checkIfExists: true)
    ]
    format='bam'

    VGAN_HAPLOCART(SAMTOOLS_BAM2FQ(input))
}
