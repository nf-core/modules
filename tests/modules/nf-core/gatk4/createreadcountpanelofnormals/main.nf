#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CREATEREADCOUNTPANELOFNORMALS } from '../../../../../modules/nf-core/gatk4/createreadcountpanelofnormals/main.nf'
include { GATK4_COLLECTREADCOUNTS             } from '../../../../../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_PREPROCESSINTERVALS           } from '../../../../../modules/nf-core/gatk4/preprocessintervals/main.nf'

workflow test_gatk4_createreadcountpanelofnormals {

    fasta = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]).collect()
    fai   = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]).collect()
    dict  = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)]).collect()

    GATK4_PREPROCESSINTERVALS ( fasta, fai, dict, [[],[]], [[],[]]).interval_list
        .map {meta,list -> list}
        .set {ch_intervals}

    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
    ],
    [
        [ id:'test2', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
    ]).combine(ch_intervals)

    GATK4_COLLECTREADCOUNTS ( input, fasta, fai, dict )

    GATK4_CREATEREADCOUNTPANELOFNORMALS (
        GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> [ [id:'test'], tsv ] })
            .groupTuple()
    )
}
