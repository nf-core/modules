#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATLAS_SPLITMERGE } from '../../../../modules/atlas/splitmerge/main.nf'
include { SAMTOOLS_INDEX as prepare_test_bai } from '../../../../modules/samtools/index/main.nf'

//Prepare the special text-file input
process prepare_test_read_group_settings {
    output:
        path 'read_group_settings.txt'
    
    script:
        """
        echo '1 paired' > read_group_settings.txt
        """
}

//prepare the tuple for ATLAS_SPLITMERGE
process prepare_test_input {
    input:
        val meta
        path bam
        path bai
        path read_group_settings
        path blacklist
    
    output:
        tuple val(meta), path(bam), path(bai), path(read_group_settings), path(blacklist)
    
    script:
        """
        """
}

//MAIN
workflow test_atlas_splitmerge {
    meta = [ id:'test', single_end:false ]
    bam = file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    bai = prepare_test_bai([
            [id:'test_bam'], 
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]).bai.map{it[1]}
    settings = prepare_test_read_group_settings()
    prepare_test_input(meta, bam, bai, settings, [])

    ATLAS_SPLITMERGE ( prepare_test_input.out )
}
