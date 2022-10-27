#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { CELLRANGER_MKGTF } from "$moduleDir/modules/nf-core/cellranger/mkgtf/main.nf"
include { CELLRANGER_MKREF } from "$moduleDir/modules/nf-core/cellranger/mkref/main.nf"
include { CELLRANGER_COUNT } from "$moduleDir/modules/nf-core/cellranger/count/main.nf"

workflow test_cellranger_count {

    input = [ [ id:'test', single_end:true, strandedness:'forward', gem: '123', samples: ["test_10x"] ], // meta map
             [ file(params.test_data['homo_sapiens']['illumina']['test_10x_1_fastq_gz'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_10x_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    CELLRANGER_MKGTF ( gtf )

    CELLRANGER_MKREF (
        fasta,
        CELLRANGER_MKGTF.out.gtf,
        reference_name
    )

    CELLRANGER_COUNT(
        input,
        CELLRANGER_MKREF.out.reference
    )
}
