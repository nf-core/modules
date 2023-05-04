#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_QC_PICARD } from '../../../../subworkflows/nf-core/bam_qc_picard/main'

workflow test_bam_qc_picard_wgs {
    input = Channel.of([ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                [],
                []
            ])
    fasta = Channel.value([
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ])
    fasta_fai = Channel.value([
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ])

    BAM_QC_PICARD ( input, fasta, fasta_fai)
}

workflow test_bam_qc_picard_targetted {
    bait        = file(params.test_data['sarscov2']['genome']['baits_interval_list'], checkIfExists: true)
    target      = file(params.test_data['sarscov2']['genome']['targets_interval_list'], checkIfExists: true)

    input = Channel.of([ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                bait,
                target
            ])
    fasta = Channel.value([
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ])
    fasta_fai = Channel.value([
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ])


    BAM_QC_PICARD ( input, fasta, fasta_fai)
}
