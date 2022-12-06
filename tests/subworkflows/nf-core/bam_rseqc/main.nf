#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_RSEQC } from '../../../../subworkflows/nf-core/bam_rseqc/main.nf'

workflow test_bam_rseqc {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    bed = file(params.test_data['sarscov2']['genome']['test_bed12'], checkIfExists: true)
    rseqc_modules = ['bam_stat', 'inner_distance', 'infer_experiment', 'junction_annotation', 'junction_saturation', 'read_distribution', 'read_duplication', 'tin']

    BAM_RSEQC ( Channel.of(input), bed, rseqc_modules )
}
