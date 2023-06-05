#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_FASTQTOBAM } from '../../../../../modules/nf-core/fgbio/fastqtobam/main.nf'

workflow test_fgbio_fastqtobam_paired_default {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    FGBIO_FASTQTOBAM ( input )
}


workflow test_fgbio_fastqtobam_paired_cram {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    FGBIO_FASTQTOBAM ( input )
}

workflow test_fgbio_fastqtobam_paired_bam {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    FGBIO_FASTQTOBAM ( input )
}

workflow test_fgbio_fastqtobam_single {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
        ]
    ]

    FGBIO_FASTQTOBAM ( input )
}

workflow test_fgbio_fastqtobam_paired_umi {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    FGBIO_FASTQTOBAM ( input )
}

workflow test_fgbio_fastqtobam_paired_custom_samplename {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    FGBIO_FASTQTOBAM ( input )
}
