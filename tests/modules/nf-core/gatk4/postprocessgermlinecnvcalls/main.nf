#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_BEDTOINTERVALLIST                                                           } from '../../../../../modules/nf-core/gatk4/bedtointervallist/main.nf'
include { GATK4_COLLECTREADCOUNTS                                                           } from '../../../../../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY as GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT } from '../../../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { GATK4_GERMLINECNVCALLER as GATK4_GERMLINECNVCALLER_COHORT                         } from '../../../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'
include { GATK4_GERMLINECNVCALLER as GATK4_GERMLINECNVCALLER_CASE                           } from '../../../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'
include { GATK4_POSTPROCESSGERMLINECNVCALLS                                                 } from '../../../../../modules/nf-core/gatk4/postprocessgermlinecnvcalls/main.nf'

workflow test_gatk4_postprocessgermlinecnvcalls {
     input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true),
    ],
    [
        [ id:'test2', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ])

    fasta = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]).collect()
    fai   = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]).collect()
    dict  = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)]).collect()

    GATK4_COLLECTREADCOUNTS ( input, fasta, fai, dict )

    bed = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true))
    contig_ploidy_table = file(params.test_data['homo_sapiens']['illumina']['contig_ploidy_priors_table'], checkIfExists: true)

    GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT (
        GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> [ [id:'test'], tsv ] })
            .groupTuple()
            .combine(bed)
            .map({ meta, counts, bed -> [ meta, counts, bed, [] ]}),
        [[],[]],
        contig_ploidy_table
    )

//
    intervals = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ]
    GATK4_BEDTOINTERVALLIST ( intervals, dict )
//
    gcnvc_cohort_input = GATK4_COLLECTREADCOUNTS.out.tsv
            .map({ meta, tsv -> [ [id:'test'], tsv ]})
            .groupTuple()
            .combine(GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT.out.calls)
            .combine(GATK4_BEDTOINTERVALLIST.out.interval_list)
            .map({ meta, counts, meta2, calls, meta3, bed -> [ meta, counts, bed, calls ]})

    GATK4_GERMLINECNVCALLER_COHORT ( gcnvc_cohort_input, [[],[]] )

    gcnvc_case_input = GATK4_COLLECTREADCOUNTS.out.tsv
            .combine(GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT.out.calls)
            .map({ meta, counts, meta2, calls -> [ [id:'test'], counts, [], calls ]})
    GATK4_GERMLINECNVCALLER_CASE ( gcnvc_case_input, GATK4_GERMLINECNVCALLER_COHORT.out.model )

    GATK4_POSTPROCESSGERMLINECNVCALLS ( GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT.out.calls,
        GATK4_GERMLINECNVCALLER_COHORT.out.model,
        GATK4_GERMLINECNVCALLER_CASE.out.calls )
}
