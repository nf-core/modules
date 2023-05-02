#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ANGSD_CONTAMINATION } from '../../../../../modules/nf-core/angsd/contamination/main.nf'
include { ANGSD_DOCOUNTS } from '../../../../../modules/nf-core/angsd/docounts/main.nf'

workflow test_angsd_contamination {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam_bai'], checkIfExists: true),
        []
    ]

    hapmap_file = [ [id:'test2'], file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/angsd/HapMapChrX.gz")]

    ANGSD_DOCOUNTS ( input )
    ANGSD_CONTAMINATION ( ANGSD_DOCOUNTS.out.icounts, hapmap_file )

    ANGSD_CONTAMINATION.out.txt.view()
}
