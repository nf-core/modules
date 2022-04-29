#!/usr/bin/env nextflow



include { SAMTOOLS_BAMTOCRAM } from '../../../../modules/samtools/bamtocram/main.nf'

workflow test_samtools_bamtocram {

    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    SAMTOOLS_BAMTOCRAM ( input, fasta, fai )
}