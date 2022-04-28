#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_QC_PICARD } from '../../../../subworkflows/nf-core/bam_qc_picard/main' addParams([:])

workflow test_bam_qc_picard_wgs {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
            ]

    BAM_QC_PICARD ( input, [], [], [] )
}

workflow test_bam_qc_picard_targetted {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]
    bait    = file(params.test_data['sarscov2']['genome']['baits_interval_list'], checkIfExists: true)
    target = file(params.test_data['sarscov2']['genome']['targets_interval_list'], checkIfExists: true)

    BAM_QC_PICARD ( input, [], bait, target )
}
