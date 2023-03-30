#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_VDJ } from '../../../../../modules/nf-core/cellranger/vdj/main.nf'

workflow test_cellranger_vdj {

    input = [ [ id:'test_cellranger_vdj_human_bcell', single_end:true, strandedness:'forward', gem: '123', samples: ["subsampled_sc5p_v2_hs_B_1k_b"] ], // meta map
             [ file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_vdj/subsampled_sc5p_v2_hs_B_1k_b_fastqs/subsampled_sc5p_v2_hs_B_1k_b_S1_L001_R1_001.fastq.gz", checkIfExists: true), 
               file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_vdj/subsampled_sc5p_v2_hs_B_1k_b_fastqs/subsampled_sc5p_v2_hs_B_1k_b_S1_L001_R2_001.fastq.gz", checkIfExists: true)
        ]
    ]

    reference_json = file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0/reference.json", checkIfExists: true)
    reference_fasta = file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0/fasta/regions.fa", checkIfExists: true)
    reference_suppfasta = file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0/fasta/supp_regions.fa", checkIfExists: true)

    reference_json.copyTo("${workDir}/vdj_reference/reference.json")
    reference_fasta.copyTo("${workDir}/vdj_reference/fasta/regions.fa")
    reference_suppfasta.copyTo("${workDir}/vdj_reference/fasta/supp_regions.fa")

    CELLRANGER_VDJ(
        input,
        file("${workDir}/vdj_reference/") // reference directory
    )
}
