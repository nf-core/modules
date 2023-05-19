#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_ARC_MKGTF } from '../../../../../modules/nf-core/cellrangerarc/mkgtf/main.nf'
include { CELLRANGER_ARC_MKREF } from '../../../../../modules/nf-core/cellrangerarc/mkref/main.nf'
include { CELLRANGER_ARC_COUNT } from '../../../../../modules/nf-core/cellrangerarc/count/main.nf'

workflow test_cellranger_arc_count {

    input = [ [ id:'test', single_end:false ], // meta map
            [  file(params.test_data['homo_sapiens']['illumina']['test_scARC_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_scARC_2_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_scARC_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_scARC_3_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_scARC_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    motifs = file(params.test_data['homo_sapiens']['genome']['genome_motif'], checkIfExists: true)
    reference_config = file(params.test_data['homo_sapiens']['genome']['genome_config'], checkIfExists: true)
    reference_name = "cellranger_arc_reference"
    lib_csv = file(params.test_data['homo_sapiens']['illumina']['multiome_lib_csv'], checkIfExists: true)

    CELLRANGER_ARC_MKGTF ( gtf )

    CELLRANGER_ARC_MKREF (
        fasta,
        CELLRANGER_MKGTF.out.gtf,
        motifs,
        reference_config,
        reference_name
    )

    CELLRANGER_ARC_COUNT(
        input,
        lib_csv,
        CELLRANGER_MKREF.out.reference
    )
}
