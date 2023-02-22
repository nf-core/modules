#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATLAS_SPLITMERGE } from '../../../../../modules/nf-core/atlas/splitmerge/main.nf'

//MAIN
workflow test_atlas_splitmerge {
    meta = [ id:'test', single_end:false ]
    bam = file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    bai = file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    settings = file(params.test_data['homo_sapiens']['illumina']['read_group_settings_txt'], checkIfExists: true)

    ATLAS_SPLITMERGE ( [meta, bam, bai, settings, []] )
}
