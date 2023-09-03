#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX as BWA_INDEX_PRIMARY} from '../../../../modules/nf-core/bwa/index/main.nf'
include { BWA_INDEX as BWA_INDEX_CONTAMINANT} from '../../../../modules/nf-core/bwa/index/main.nf'
include { FASTQ_ALIGN_BAMCMP_BWA } from '../../../../subworkflows/nf-core/fastq_align_bamcmp_bwa/main.nf'

workflow test_fastq_align_bamcmp_bwa {

    input = [
        [ id:'test'], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]

    primary_fasta = [
        [ id:'homo_sapiens_genome'], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    contaminant_fasta = [
        [ id:'sarscov2_genome'], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BWA_INDEX_PRIMARY ( primary_fasta )
    BWA_INDEX_CONTAMINANT ( contaminant_fasta )
    FASTQ_ALIGN_BAMCMP_BWA ( input, BWA_INDEX_PRIMARY.out.index, BWA_INDEX_CONTAMINANT.out.index, false, primary_fasta )

}
