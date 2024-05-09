#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HYPO } from '../../../../modules/nf-core/hypo/main.nf'

workflow test_hypo {
    
    bam_sr = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    reads = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    draft = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    genome_size = 29903
    reads_coverage = 5
    HYPO ( bam_sr, reads, draft, genome_size, reads_coverage )
}
