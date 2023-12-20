#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGERARC_MKREF } from '../../../../../modules/nf-core/cellrangerarc/mkref/main.nf'
include { CELLRANGERARC_COUNT } from '../../../../../modules/nf-core/cellrangerarc/count/main.nf'
include { UNZIP } from '../../../../../modules/nf-core/unzip/main.nf'

workflow test_cellrangerarc_count {

    input = [ [ id:'test' ], [ 'gex','atac' ], [ 'SRR18907480_chr19_sub','SRR18907481_chr19_sub' ],
            [  file(params.test_data['mus_musculus']['illumina']['test_scARC_gex_R1_fastq_gz'], checkIfExists: true),
                file(params.test_data['mus_musculus']['illumina']['test_scARC_gex_R2_fastq_gz'], checkIfExists: true),
                file(params.test_data['mus_musculus']['illumina']['test_scARC_atac_R1_fastq_gz'], checkIfExists: true),
                file(params.test_data['mus_musculus']['illumina']['test_scARC_atac_R2_fastq_gz'], checkIfExists: true),
                file(params.test_data['mus_musculus']['illumina']['test_scARC_atac_I2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = [ [], file(params.test_data['mus_musculus']['genome']['genome_19_fasta'], checkIfExists: true)]
    gtf = file(params.test_data['mus_musculus']['genome']['genome_19_gtf'], checkIfExists: true)
    motifs = file(params.test_data['homo_sapiens']['genome']['genome_motifs'], checkIfExists: true)
    reference_config = file(params.test_data['mus_musculus']['illumina']['genome_config'], checkIfExists: true)
    reference_name = "cellrangerarc_reference"
    lib_csv = file(params.test_data['mus_musculus']['illumina']['multiome_lib_csv'], checkIfExists: true)

    UNZIP( fasta )

    CELLRANGERARC_MKREF ( UNZIP.out.unzipped_archive.map{ it[1] } + "/chr19.fa",
                            gtf,
                            motifs,
                            reference_config,
                            reference_name )

    CELLRANGERARC_COUNT(
        input,
        CELLRANGERARC_MKREF.out.reference
    )
}
