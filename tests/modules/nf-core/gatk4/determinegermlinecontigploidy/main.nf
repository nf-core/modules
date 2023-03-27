#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_DETERMINEGERMLINECONTIGPLOIDY as GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT } from '../../../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY as GATK4_DETERMINEGERMLINECONTIGPLOIDY_CASE   } from '../../../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { GATK4_COLLECTREADCOUNTS                                                           } from '../../../../../modules/nf-core/gatk4/collectreadcounts/main.nf'

workflow test_gatk4_determinegermlinecontigploidy_case {

    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ],
    [
        [ id:'test2', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ])

    fasta = []
    fai = []
    dict = []

    GATK4_COLLECTREADCOUNTS (
        input,
        fasta,
        fai,
        dict
    )

    bed = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true))
    contig_ploidy_table = file(params.test_data['homo_sapiens']['illumina']['contig_ploidy_priors_table'], checkIfExists: true)

    GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT (
        GATK4_COLLECTREADCOUNTS.out.tsv.tap { case_input }
            .map({ meta, tsv -> [ [id:'test'], tsv ] })
            .groupTuple()
            .combine(bed)
            .map({ meta, counts, bed -> [ meta, counts, bed, [] ]}),
        contig_ploidy_table,
        []
    )

    GATK4_DETERMINEGERMLINECONTIGPLOIDY_CASE (
        case_input
            .first()
            .map({ meta, tsv -> [ [id:'test_case'], tsv, [], [] ] }),
        [],
        GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT.out.model
            .map({ meta, model -> model})
    )
}

workflow test_gatk4_determinegermlinecontigploidy_cohort {

    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ],
    [
        [ id:'test2', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ])

    fasta = []
    fai = []
    dict = []

    GATK4_COLLECTREADCOUNTS (
        input,
        fasta,
        fai,
        dict
    )

    contig_ploidy_table = file(params.test_data['homo_sapiens']['illumina']['contig_ploidy_priors_table'], checkIfExists: true)

    GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT (
        GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> [ [id:'test'], tsv ] })
            .groupTuple()
            .map({ meta, counts -> [ meta, counts, [], [] ]}),
        contig_ploidy_table,
        []
    )
}

workflow test_gatk4_determinegermlinecontigploidy_cohort_include_intervals {

    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ],
    [
        [ id:'test2', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ])

    fasta = []
    fai = []
    dict = []

    GATK4_COLLECTREADCOUNTS (
        input,
        fasta,
        fai,
        dict
    )

    bed = Channel.of(file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true))
    contig_ploidy_table = file(params.test_data['homo_sapiens']['illumina']['contig_ploidy_priors_table'], checkIfExists: true)

    GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT (
        GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> [ [id:'test'], tsv ] })
            .groupTuple()
            .combine(bed)
            .map({ meta, counts, bed -> [ meta, counts, bed, [] ]}),
        contig_ploidy_table,
        []
    )
}

workflow test_gatk4_determinegermlinecontigploidy_cohort_exclude_intervals {

    bed_1 = Channel.of(file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true))
            .map({bed -> return bed.text.replace("40001","10001") + "chr22\t10001\t40001\n" })
            .collectFile( name: "genome.bed" )

    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ],
    [
        [ id:'test2', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
    ]).combine(bed_1)

    fasta = []
    fai = []
    dict = []

    GATK4_COLLECTREADCOUNTS (
        input,
        fasta,
        fai,
        dict
    )

    bed_2 = Channel.of(file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true))
            .map({bed -> return bed.text.replaceFirst("0","10001") })
            .collectFile( name: "genome.bed" )
    contig_ploidy_table = file(params.test_data['homo_sapiens']['illumina']['contig_ploidy_priors_table'], checkIfExists: true)

    GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT (
        GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> [ [id:'test'], tsv ] })
            .groupTuple()
            .combine(bed_2)
            .map({ meta, counts, bed -> [ meta, counts, [], bed ]}),
        contig_ploidy_table,
        []
    )
}
