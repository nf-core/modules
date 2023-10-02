#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_DOCOUNTS_CONTAMINATION_ANGSD } from '../../../../subworkflows/nf-core/bam_docounts_contamination_angsd/main.nf'

workflow test_bam_docounts_contamination_angsd {

    input = [
        [
            [ id:'test', single_end: false ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'], checkIfExists: true)
        ]
    ]

    bai = [
        [
            [ id:'test', single_end: false ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam_bai'], checkIfExists: true)
        ]
    ]

    hapmap_file = [
        [
            [ id:'test2' ], // meta map
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/angsd/HapMapChrX.gz", checkIfExists: true)
        ]
    ]

    BAM_DOCOUNTS_CONTAMINATION_ANGSD ( Channel.fromList(input), Channel.fromList(bai), Channel.fromList(hapmap_file) )
}
