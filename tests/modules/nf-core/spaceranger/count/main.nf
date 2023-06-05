#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPACERANGER_MKGTF } from '../../../../../modules/nf-core/spaceranger/mkgtf/main.nf'
include { SPACERANGER_MKREF } from '../../../../../modules/nf-core/spaceranger/mkref/main.nf'
include {
    SPACERANGER_COUNT as SPACERANGER_COUNT_FFPE_V1;
    SPACERANGER_COUNT as SPACERANGER_COUNT_FFPE_CYTASSIST

} from '../../../../../modules/nf-core/spaceranger/count/main.nf'

workflow test_spaceranger_count {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    SPACERANGER_MKGTF ( gtf )

    SPACERANGER_MKREF (
        fasta,
        SPACERANGER_MKGTF.out.gtf,
        reference_name
    )

    ch_spaceranger_ref = SPACERANGER_MKREF.out.reference

    SPACERANGER_COUNT_FFPE_V1 (
        [
            [
                id: "Visium_FFPE_Human_Ovarian_Cancer",
                slide: "V10L13-020",
                area: "D1"
            ],
            [
                file(params.test_data['homo_sapiens']['10xgenomics']['test_10x_ffpe_v1_fastq_1_gz']),
                file(params.test_data['homo_sapiens']['10xgenomics']['test_10x_ffpe_v1_fastq_2_gz'])
            ],
            file(params.test_data['homo_sapiens']['10xgenomics']['test_10x_ffpe_v1_image']),
            [],
            file('https://s3.us-west-2.amazonaws.com/10x.spatial-slides/gpr/V10L13/V10L13-020.gpr')
        ],
        ch_spaceranger_ref,
        [],
    )

    SPACERANGER_COUNT_FFPE_CYTASSIST (
        [
            [
                id: "CytAssist_11mm_FFPE_Human_Glioblastoma_2",
                slide: "V52Y10-317",
                area: "B1"
            ],
            [
                file(params.test_data['homo_sapiens']['10xgenomics']['test_10x_ffpe_cytassist_fastq_1_gz']),
                file(params.test_data['homo_sapiens']['10xgenomics']['test_10x_ffpe_cytassist_fastq_2_gz'])
            ],
            file(params.test_data['homo_sapiens']['10xgenomics']['test_10x_ffpe_cytassist_image']),
            [],
            file('https://s3.us-west-2.amazonaws.com/10x.spatial-slides/gpr/V52Y10/V52Y10-317.gpr')
        ],
        ch_spaceranger_ref,
        file(params.test_data['homo_sapiens']['10xgenomics']['test_10x_ffpe_cytassist_probeset']),
    )
}
