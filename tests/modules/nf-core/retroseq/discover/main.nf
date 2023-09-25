#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RETROSEQ_DISCOVER } from '../../../../../modules/nf-core/retroseq/discover/main.nf'

workflow test_retroseq_discover {

    te_elements = new File("elements.txt")
    te_elements.write("ALU\tannotated.bed")

    bam = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
    ]

    extra_files = [
        file(params.test_data['homo_sapiens']['genome']['genome_21_annotated_bed'], checkIfExists: true)
    ]

    RETROSEQ_DISCOVER(bam, file(te_elements), extra_files)
}
