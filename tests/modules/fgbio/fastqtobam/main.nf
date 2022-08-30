#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_FASTQTOBAM } from '../../../../modules/fgbio/fastqtobam/main.nf'
include { FGBIO_FASTQTOBAM as FGBIO_FASTQTOBAM_UMI  } from '../../../../modules/fgbio/fastqtobam/main.nf'
include { FGBIO_FASTQTOBAM as FGBIO_FASTQTOBAM_CUSTOM_SAMPLENAME } from '../../../../modules/fgbio/fastqtobam/main.nf'

workflow test_fgbio_fastqtobam_paired {

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
    read_structure = "+T 12M11S+T"

    FGBIO_FASTQTOBAM_UMI ( input )
}

workflow test_fgbio_fastqtobam_paired_custom_samplename {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    read_structure = "+T 12M11S+T"

    FGBIO_FASTQTOBAM_CUSTOM_SAMPLENAME ( input )
}

