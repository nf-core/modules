#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_STATS } from '../../../../modules/samtools/stats/main.nf' addParams( options: [:] )

workflow test_samtools_stats {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]

    SAMTOOLS_STATS ( input, [])
}

workflow test_samtools_stats_cram {
   input = [ [ id: 'test', single_end:true ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true)
            ]
    fasta   = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    SAMTOOLS_STATS ( input, fasta )
}
