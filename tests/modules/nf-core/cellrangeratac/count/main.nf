#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGERATAC_MKREF } from '../../../../../modules/nf-core/cellrangeratac/mkref/main.nf'
include { CELLRANGERATAC_COUNT } from '../../../../../modules/nf-core/cellrangeratac/count/main.nf'

workflow test_cellrangeratac_count {

    input = [ [ id:'test', single_end:false, samples: ["test_scATAC"] ], // meta map
             //[  file(params.test_data['homo_sapiens']['illumina']['test_scATAC_1_fastq_gz'], checkIfExists: true),
             //   file(params.test_data['homo_sapiens']['illumina']['test_scATAC_3_fastq_gz'], checkIfExists: true),
             //   file(params.test_data['homo_sapiens']['illumina']['test_scATAC_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_scATAC_gtf'], checkIfExists: true)
    motifs = file(params.test_data['homo_sapiens']['genome']['genome_motif'], checkIfExists: true)
    reference_config = file(params.test_data['homo_sapiens']['genome']['genome_config'], checkIfExists: true)
    reference_name = "cellrangeratac_reference"

    CELLRANGERATAC_MKREF (
        fasta,
        gtf,
        motifs,
        reference_config,
        reference_name
    )

    CELLRANGERATAC_COUNT(
        input,
        CELLRANGER_MKREF.out.reference
    )
}
