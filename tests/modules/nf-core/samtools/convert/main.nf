#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_CONVERT as SAMTOOLS_BAMTOCRAM } from '../../../../../modules/nf-core/samtools/convert/main.nf'
include { SAMTOOLS_CONVERT as SAMTOOLS_CRAMTOBAM } from '../../../../../modules/nf-core/samtools/convert/main.nf'

workflow test_samtools_convert_bamtocram {

    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    SAMTOOLS_BAMTOCRAM ( input, fasta, fai )
}

workflow test_samtools_convert_cramtobam {

    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true)
            ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    SAMTOOLS_CRAMTOBAM ( input, fasta, fai )
}
