#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPVARIANT } from '../../../modules/deepvariant/main.nf' addParams( options: ["args": "--regions=\"chr20:10,000,000-10,010,000\" --model_type=WGS"] )

workflow test_deepvariant {

    bam_tuple_ch = Channel.of([[ id:'test', single_end:false ], // meta map
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    DEEPVARIANT ( bam_tuple_ch, fasta, fai)
}
