#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_VIEW } from '../../../../modules/samtools/view/main.nf'

workflow test_samtools_view {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]

    SAMTOOLS_VIEW ( input, [] )
}

workflow test_samtools_view_cram {
   input = [ [ id: 'test' ], // meta map
               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true)
            ]
    fasta   = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    SAMTOOLS_VIEW ( input, fasta )
}
