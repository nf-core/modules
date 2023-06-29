#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGERATAC_MKREF } from '../../../../../modules/nf-core/cellrangeratac/mkref/main.nf'
include { CELLRANGERATAC_COUNT } from '../../../../../modules/nf-core/cellrangeratac/count/main.nf'
include { UNZIP } from '../../../../../modules/nf-core/unzip/main.nf'

workflow test_cellrangeratac_count {

    input = [ [ id:'test', single_end:false, samples: ["test_scATAC"] ], // meta map
                 [  file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_scATAC_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_scATAC_3_fastq_gz'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_scATAC_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = [ [], file(params.test_data['homo_sapiens']['genome']['genome_1_fasta'], checkIfExists: true) ]
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_1_gtf'], checkIfExists: true)
    motifs = file(params.test_data['homo_sapiens']['genome']['genome_motifs'], checkIfExists: true)
    reference_config = file(params.test_data['homo_sapiens']['genome']['genome_config'], checkIfExists: true)
    reference_name = "cellrangeratac_reference"

    UNZIP( fasta )

    CELLRANGERATAC_MKREF (
        UNZIP.out.unzipped_archive.map { it[1] } + "/genome.fasta",
        gtf,
        motifs,
        reference_config,
        reference_name
    )

    CELLRANGERATAC_COUNT(
        input,
        CELLRANGERATAC_MKREF.out.reference
    )
}
