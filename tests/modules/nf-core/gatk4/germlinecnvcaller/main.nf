#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_COLLECTREADCOUNTS } from '../../../../../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY as GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT; GATK4_DETERMINEGERMLINECONTIGPLOIDY as GATK4_DETERMINEGERMLINECONTIGPLOIDY_CASE } from '../../../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { UNTAR as UNTAR_PLOIDY_COHORT; UNTAR as UNTAR_PLOIDY_CASE; UNTAR as UNTAR_MODEL } from '../../../../../modules/nf-core/untar/main.nf'
include { GATK4_GERMLINECNVCALLER as GATK4_GERMLINECNVCALLER_COHORT; GATK4_GERMLINECNVCALLER as GATK4_GERMLINECNVCALLER_CASE } from '../../../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'

workflow test_gatk4_germlinecnvcaller_case {
    intervals = file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    bed = Channel.of(file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true))
    priors =  file(params.test_data['homo_sapiens']['illumina']['contig_ploidy_priors_table'], checkIfExists: true)
    inputs = Channel.of([
        [ id:'test1' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        intervals
    ],
    [
        [ id:'test2' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        intervals
    ])
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_COLLECTREADCOUNTS ( inputs, fasta, fai, dict )
    dgcp_cohort_input = GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> [ [id:'test'], tsv ] })
            .groupTuple()
            .combine(bed)
            .map({ meta, counts, bed -> [ meta, counts, bed, [] ]})
    input = GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> return [[id:"test"], tsv ]})
            .groupTuple()

    GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT ( dgcp_cohort_input, priors, [] )

    model = GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT.out.model
            .map({ meta, model -> model})
    dgcp_case_input = dgcp_cohort_input
            .map({ meta, counts, bed, exclude_beds -> [ meta, counts, [], [] ]})

    GATK4_DETERMINEGERMLINECONTIGPLOIDY_CASE ( dgcp_case_input, [], model )
    GATK4_GERMLINECNVCALLER_COHORT ( input, intervals, [], GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT.out.calls.collect{ it[1] } )
    GATK4_GERMLINECNVCALLER_CASE ( input, [], GATK4_GERMLINECNVCALLER_COHORT.out.model.collect{ it[1] }, GATK4_DETERMINEGERMLINECONTIGPLOIDY_CASE.out.calls.collect{ it[1] } )
}

workflow test_gatk4_germlinecnvcaller_cohort {
    intervals = file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    bed = Channel.of(file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true))
    priors =  file(params.test_data['homo_sapiens']['illumina']['contig_ploidy_priors_table'], checkIfExists: true)
    inputs = Channel.of([
        [ id:'test1' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        intervals
    ],
    [
        [ id:'test2' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        intervals
    ])
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_COLLECTREADCOUNTS ( inputs, fasta, fai, dict )
    dgcp_cohort_input = GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> [ [id:'test'], tsv ] })
            .groupTuple()
            .combine(bed)
            .map({ meta, counts, bed -> [ meta, counts, bed, [] ]})
    input = GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> return [[id:"test"], tsv ]})
            .groupTuple()
    model = []

    GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT ( dgcp_cohort_input, priors, [] )
    GATK4_GERMLINECNVCALLER_COHORT ( input, intervals, model, GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT.out.calls.collect{ it[1] } )
}
