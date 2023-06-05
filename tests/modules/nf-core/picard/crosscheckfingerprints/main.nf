#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_CROSSCHECKFINGERPRINTS } from '../../../../../modules/nf-core/picard/crosscheckfingerprints/main.nf'

workflow test_picard_crosscheckfingerprints {

    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)],
    ]
    PICARD_CROSSCHECKFINGERPRINTS ( input,[], file(params.test_data['homo_sapiens']['genome']['haplotype_map'], checkIfExists: true))
}
