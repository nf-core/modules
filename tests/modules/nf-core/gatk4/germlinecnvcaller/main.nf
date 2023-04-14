#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_COLLECTREADCOUNTS                                                           } from '../../../../../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY as GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT } from '../../../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY as GATK4_DETERMINEGERMLINECONTIGPLOIDY_CASE   } from '../../../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { GATK4_GERMLINECNVCALLER as GATK4_GERMLINECNVCALLER_COHORT                         } from '../../../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'
include { GATK4_GERMLINECNVCALLER as GATK4_GERMLINECNVCALLER_CASE                           } from '../../../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'

workflow test_gatk4_germlinecnvcaller_case {
    bed = Channel.of(file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true))
    priors =  file(params.test_data['homo_sapiens']['illumina']['contig_ploidy_priors_table'], checkIfExists: true)
    inputs = Channel.of([
        [ id:'test1' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ],
    [
        [ id:'test2' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ])
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    ploidy_calls = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_ploidy_calls'], checkIfExists: true) ])
    cnv_model = Channel.of(file(params.test_data['homo_sapiens']['genome']['genome_germline_cnv_model'], checkIfExists: true))

    GATK4_COLLECTREADCOUNTS ( inputs, fasta, fai, dict )

    gcnvc_case_input = GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> [ [id:'test'], tsv ]})
            .groupTuple()
            .map({ meta, counts -> [ meta, counts, [] ]})
    GATK4_GERMLINECNVCALLER_CASE ( gcnvc_case_input, cnv_model, ploidy_calls.collect{ it[1] } )
    GATK4_GERMLINECNVCALLER_CASE.out.calls
}

workflow test_gatk4_germlinecnvcaller_cohort {
    bed = Channel.of(file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true))
    priors =  file(params.test_data['homo_sapiens']['illumina']['contig_ploidy_priors_table'], checkIfExists: true)
    inputs = Channel.of([
        [ id:'test1' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ],
    [
        [ id:'test2' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ])
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    ploidy_calls = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_ploidy_calls'], checkIfExists: true) ])

    GATK4_COLLECTREADCOUNTS ( inputs, fasta, fai, dict )

    gcnvc_cohort_input = GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> return [[id:"test"], tsv ]})
            .groupTuple()
            .combine(bed)
            .map({ meta, counts, bed -> [ meta, counts, bed ]})
    GATK4_GERMLINECNVCALLER_COHORT ( gcnvc_cohort_input, [], ploidy_calls.collect{ it[1] } )
    GATK4_GERMLINECNVCALLER_COHORT.out.calls
    GATK4_GERMLINECNVCALLER_COHORT.out.model
}
