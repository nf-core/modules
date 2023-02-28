#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT_BATCH as CNVKIT_HYBRID    } from '../../../../../modules/nf-core/cnvkit/batch/main.nf'
include { CNVKIT_BATCH as CNVKIT_WGS       } from '../../../../../modules/nf-core/cnvkit/batch/main.nf'
include { CNVKIT_BATCH as CNVKIT_TUMORONLY } from '../../../../../modules/nf-core/cnvkit/batch/main.nf'
include { CNVKIT_BATCH as CNVKIT_GERMLINE  } from '../../../../../modules/nf-core/cnvkit/batch/main.nf'
include { CNVKIT_BATCH as CNVKIT_PON       } from '../../../../../modules/nf-core/cnvkit/batch/main.nf'

workflow test_cnvkit_hybrid_somatic {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true)
    ]
    fasta   = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    targets = file(params.test_data['sarscov2']['genome']['baits_bed'], checkIfExists: true)

    CNVKIT_HYBRID ( input, fasta, [], targets, [], false )
}

workflow test_cnvkit_wgs_somatic {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    CNVKIT_WGS ( input, fasta, [], [], [], false )
}

workflow test_cnvkit_cram_wgs_somatic {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true)
    ]
    fasta     = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai     = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    CNVKIT_WGS ( input, fasta, fasta_fai, [], [], false )
}


workflow test_cnvkit_tumoronly_hybrid_bam {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        []
    ]
    reference = file(params.test_data['homo_sapiens']['genome']['genome_21_reference_cnn'], checkIfExists: true)

    CNVKIT_TUMORONLY ( input, [], [],  [], reference, false )
}

workflow test_cnvkit_tumoronly_hybrid_cram {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        []
    ]
    fasta     = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    reference = file(params.test_data['homo_sapiens']['genome']['genome_21_reference_cnn'], checkIfExists: true)

    CNVKIT_TUMORONLY ( input, fasta, [], [], reference, false )
}

workflow test_cnvkit_germline_hybrid_cram {

    input = [
        [ id:'test'], // meta map
        [],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'], checkIfExists: true)
    ]
    fasta     = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    targets   = file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)

    CNVKIT_GERMLINE ( input, fasta, fasta_fai, targets,  [], false )
}

workflow test_cnvkit_germline_hybrid_bam {

    input = [
        [ id:'test'], // meta map
        [],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true)
    ]
    fasta     = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    targets   = file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)

    CNVKIT_GERMLINE ( input, fasta, [], targets,  [], false )
}

workflow test_cnvkit_pon {

    input = [
        [ id:'test'], // meta map
        [],
        [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    CNVKIT_PON ( input, fasta, [], [], [], true )
}
