#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_DOCOUNTS_CONTAMINATION_ANGSD } from '../../../../subworkflows/nf-core/bam_docounts_contamination_angsd/main.nf'

workflow test_bam_docounts_contamination_angsd {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam_bai'], checkIfExists: true),
        hapmap_file = file("https://github.com/jbv2/nf-core-test-datasets/raw/modules/data/delete_me/angsd/HapMapChrX.gz"), checkIfExists: true
    ]

    BAM_DOCOUNTS_CONTAMINATION_ANGSD ( input )
}
