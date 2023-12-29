#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GRIDSS_GRIDSS } from '../../../../../modules/nf-core/gridss/gridss/main.nf'
include { BWA_INDEX     } from '../../../../../modules/nf-core/bwa/index/main.nf'
include { GRIDSS_GRIDSSGENERATEPONBEDPE } from '../../../../../modules/nf-core/gridss/gridssgenerateponbedpe/main.nf'


workflow test_gridss_gridssgenerateponbedpe {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        []
    ]
    fasta = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[id:'fasta_fai'],file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]

    bwa_index = BWA_INDEX(fasta).index

    GRIDSS_GRIDSS(input,fasta,fasta_fai,bwa_index)
    input2  = GRIDSS_GRIDSS.out.vcf.map{ it -> tuple( it[0], it[1], [],[] )}
    GRIDSS_GRIDSSGENERATEPONBEDPE ( input2,fasta,fasta_fai,bwa_index )
}
