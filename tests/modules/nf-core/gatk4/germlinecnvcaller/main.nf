#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_COLLECTREADCOUNTS } from '../../../../../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_COLLECTREADCOUNTS as GATK4_COLLECTREADCOUNTS_INPUT1 } from '../../../../../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_COLLECTREADCOUNTS as GATK4_COLLECTREADCOUNTS_INPUT2 } from '../../../../../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../../../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { UNTAR } from '../../../../../modules/nf-core/untar/main.nf'
include { GATK4_GERMLINECNVCALLER } from '../../../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'

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
    model = []
    
    
    GATK4_COLLECTREADCOUNTS ( inputs, fasta, fai, dict )
    dgcp_input = GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> [ [id:'test'], tsv ] })
            .groupTuple()
            .combine(bed)
            .map({ meta, counts, bed -> [ meta, counts, bed, [] ]})
    read_counts = GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv  -> return [[id:"test"], tsv ]})
            .groupTuple()
    GATK4_DETERMINEGERMLINECONTIGPLOIDY ( dgcp_input, priors, [] )
    UNTAR (GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls)
    GATK4_GERMLINECNVCALLER ( read_counts, intervals, model, UNTAR.out.untar.collect{ it[1] } )
}
