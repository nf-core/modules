#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNIVERSC_MKGTF } from '../../../../../modules/nf-core/universc/cellranger_os_mkgtf/main.nf'
include { UNIVERSC_MKREF } from '../../../../../modules/nf-core/universc/cellranger_os_mkref/main.nf'
include { UNIVERSC_COUNT } from '../../../../../modules/nf-core/universc/cellranger_os_count/main.nf'

workflow test_cellranger_count {

    input = [ [ id:'123', chemistry:'SC3Pv3', single_end:false, strandedness:'forward', gem: '123', samples: ["test_10x"] ], // meta map
             [ file(params.test_data['homo_sapiens']['illumina']['test_10x_1_fastq_gz'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_10x_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    UNIVERSC_MKGTF ( gtf )

    UNIVERSC_MKREF (
        fasta,
        UNIVERSC_MKGTF.out.gtf,
        reference_name
    )

    UNIVERSC_COUNT (
        input,
        UNIVERSC_MKREF.out.reference
    )
}
