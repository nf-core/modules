#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_VDJ } from '../../../../../modules/nf-core/cellranger/vdj/main.nf'

workflow test_cellranger_vdj {

    input = [ [ id:'subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ', single_end:false, strandedness:'auto' ], // meta map
        [ file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_vdj_fastq_1_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_vdj_fastq_2_gz'], checkIfExists: true)
        ],
    ]

    reference_json = file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_vdj_ref_json'], checkIfExists: true)
    reference_fasta = file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_vdj_ref_fasta'], checkIfExists: true)
    reference_suppfasta = file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_vdj_ref_suppfasta'], checkIfExists: true)

    reference_json.copyTo("${workDir}/vdj_reference/reference.json")
    reference_fasta.copyTo("${workDir}/vdj_reference/fasta/regions.fa")
    reference_suppfasta.copyTo("${workDir}/vdj_reference/fasta/supp_regions.fa")

    CELLRANGER_VDJ(
        input,
        file("${workDir}/vdj_reference/") // reference directory
    )
}
