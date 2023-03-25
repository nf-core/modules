#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_NGSCHECKMATE } from '../../../../subworkflows/nf-core/bam_ngscheckmate/main.nf'
include { BEDTOOLS_MAKEWINDOWS } from '../../../../modules/nf-core/bedtools/makewindows/main.nf'

workflow test_bam_ngscheckmate {

    inputBed = [ [ id:'test_bed'],
                file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)]

    BEDTOOLS_MAKEWINDOWS(Channel.of(inputBed))
    snp_bed  = BEDTOOLS_MAKEWINDOWS.out.bed

    input   = [[
            [ id:'test1' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
        ],[
            [ id:'test2' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam'], checkIfExists: true)
        ]
    ]

    fasta    = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    ch_input = Channel.fromList(input)
    ch_snp_bed = BEDTOOLS_MAKEWINDOWS.out.bed.map{it[1]}
    ch_fasta = Channel.of(fasta)

    BAM_NGSCHECKMATE( ch_input, ch_snp_bed , ch_fasta)

}
