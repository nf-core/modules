#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_BEDTOINTERVALLIST } from '../../../../modules/gatk4/bedtointervallist/main.nf'
include { GATK4_INTERVALLISTTOOLS } from '../../../../modules/gatk4/intervallisttools/main.nf'

// include { GATK4_PREPROCESSINTERVALS } from '../../../../modules/gatk4/preprocessintervals/main.nf'
// include { GATK4_ANNOTATEINTERVALS } from '../../../../modules/gatk4/annotateintervals/main.nf'
// include { GATK4_COLLECTREADCOUNTS } from '../../../../modules/gatk4/collectreadcounts/main.nf'
// include { GATK4_FILTERINTERVALS } from '../../../../modules/gatk4/filterintervals/main.nf'

workflow test_gatk4_intervallisttools {

    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_BEDTOINTERVALLIST ( input, dict )
    GATK4_INTERVALLISTTOOLS ( GATK4_BEDTOINTERVALLIST.out.interval_list )
}

// workflow test_gatk4_intervallisttools {
    
//     bed = file(params.test_data['homo_sapiens']['genome']['genome_blacklist_interval_bed'], checkIfExists: true)
//     fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
//     fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
//     dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
//     input = [
//         [ id:'test' ], // meta map
//         file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
//         file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
//     ]

//     GATK4_PREPROCESSINTERVALS ( bed, fasta, fai, dict )
//     GATK4_ANNOTATEINTERVALS ( GATK4_PREPROCESSINTERVALS.out.interval_list, fasta, fai, dict  )
//     GATK4_COLLECTREADCOUNTS ( input, GATK4_PREPROCESSINTERVALS.out.interval_list, fasta, fai, dict )
//     GATK4_FILTERINTERVALS ( GATK4_COLLECTREADCOUNTS.out.tsv, GATK4_PREPROCESSINTERVALS.out.interval_list, GATK4_ANNOTATEINTERVALS.out.tsv )
//     GATK4_INTERVALLISTTOOLS ( GATK4_FILTERINTERVALS.out.interval_list )
// }
