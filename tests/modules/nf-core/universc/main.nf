#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKGTF } from '../../../../modules/nf-core/cellranger/mkgtf/main.nf'
include { CELLRANGER_MKREF } from '../../../../modules/nf-core/cellranger/mkref/main.nf'
include { UNIVERSC } from '../../../../modules/nf-core/universc/main.nf'

workflow test_universc_10x {

    input = [
        [
            id: '123',
            technology: '10x',
            chemistry: 'SC3Pv3',
            single_end: false,
            strandedness: 'forward',
            samples: ["test_10x"]
        ], // meta map
        [
            file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_gex1_fastq_1_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_gex1_fastq_2_gz'], checkIfExists: true)
        ]
    ]

    fasta          = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf            = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    CELLRANGER_MKGTF ( gtf )

    CELLRANGER_MKREF (
        fasta,
        CELLRANGER_MKGTF.out.gtf,
        reference_name
    )

    UNIVERSC (
        input,
        CELLRANGER_MKREF.out.reference
    )
}
