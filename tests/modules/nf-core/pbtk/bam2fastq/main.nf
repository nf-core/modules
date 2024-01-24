#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PBTK_BAM2FASTQ } from '../../../../../modules/nf-core/pbtk/bam2fastq/main.nf'

workflow test_pbtk_bam2fastq {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.modules_testdata_base_path + 'genomics/homo_sapiens/pacbio/bam/alz.bam', checkIfExists: true),
        file(params.modules_testdata_base_path + 'genomics/homo_sapiens/pacbio/bam/alz.bam.pbi', checkIfExists: true),
    ]

    PBTK_BAM2FASTQ ( input )
}
