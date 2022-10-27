#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_ATAC_MKREF } from '../../../../../modules/nf-core/cellranger/mkref/main.nf'
include { CELLRANGER_ATAC_COUNT } from '../../../../../modules/nf-core/cellranger/count/main.nf'

workflow test_cellranger_atac_count {

    input = [ [ id:'test', single_end:true, strandedness:'forward', gem: '123', samples: ["test_10x"] ], // meta map
             [  file(params.test_data['homo_sapiens']['illumina']['test_10x_scATAC_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_10x_scATAC_3_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_10x_scATAC_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    motif = file(params.test_data['homo_sapiens']['genome']['genome_motif'], checkIfExists: true)
    reference_config = file(, checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    CELLRANGER_ATAC_MKREF (
        fasta,
        gtf,
        motif,
        reference_config,
        reference_name
    )

    CELLRANGER_ATAC_COUNT(
        input,
        CELLRANGER_MKREF.out.reference
    )
}
