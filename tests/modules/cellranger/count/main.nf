#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKREF } from '../../../../modules/cellranger/mkref/main.nf' addParams( options: [:] )
include { CELLRANGER_COUNT } from '../../../../modules/cellranger/count/main.nf' addParams( options: [ args: '--chemistry SC3Pv3'] )

workflow test_cellranger_count {

    input = [ [ id:'test', single_end:true, strandedness:'forward', gem: '123', samples: ["test_10x"] ], // meta map
             [ file(params.test_data['homo_sapiens']['illumina']['test_10x_1_fastq_gz'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_10x_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    reference_name = "homo_sapiens_chr22_reference"

    CELLRANGER_MKREF (
        fasta,
        gtf,
        reference_name
    )

    CELLRANGER_COUNT(
        input,
        CELLRANGER_MKREF.out.reference
    )
}
