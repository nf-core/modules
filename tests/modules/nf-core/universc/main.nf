#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKGTF_OS } from '../../../../modules/nf-core/universc/mkgtf/main.nf'
include { CELLRANGER_MKREF_OS } from '../../../../modules/nf-core/universc/mkref/main.nf'
include { CELLRANGER_COUNT_OS } from '../../../../modules/nf-core/universc/main.nf'
include { UNIVERSC } from '../../../../modules/nf-core/universc/main.nf'

workflow test_cellranger_10x {
    
    input = [ [ id:'test', chemistry:'SC3Pv3', single_end:false, strandedness:'forward', gem: '123', samples: ["test_10x"] ], // meta map
             [ file(params.test_data['homo_sapiens']['illumina']['test_10x_1_fastq_gz'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_10x_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    CELLRANGER_MKGTF_OS ( gtf )

    CELLRANGER_MKREF_OS (
        fasta,
        CELLRANGER_MKGTF_OS.out.gtf,
        reference_name
    )

    CELLRANGER_COUNT_OS (
        input,
        CELLRANGER_MKREF_OS.out.reference
    )
}

workflow test_universc_10x {
    
    input = [ [ id:'123', technology:'10x', chemistry:'SC3Pv3', single_end:false, strandedness:'forward', samples: ["test_10x"] ], // meta map
             [ file(params.test_data['homo_sapiens']['illumina']['test_10x_1_fastq_gz'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_10x_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    CELLRANGER_MKGTF_OS ( gtf )

    CELLRANGER_MKREF_OS (
        fasta,
        CELLRANGER_MKGTF_OS.out.gtf,
        reference_name
    )

    UNIVERSC (
        input,
        CELLRANGER_MKREF_OS.out.reference
    )
}
