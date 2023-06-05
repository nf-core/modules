#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATLAS_RECAL } from '../../../../../modules/nf-core/atlas/recal/main.nf'

workflow test_atlas_recal {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        [],
        []
    ]
    alleles = []
    invariant_sites = []

    ATLAS_RECAL ( input, alleles, invariant_sites )
}


