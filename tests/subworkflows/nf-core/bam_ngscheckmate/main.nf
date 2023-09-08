#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_NGSCHECKMATE } from '../../../../subworkflows/nf-core/bam_ngscheckmate/main.nf'
include { BEDTOOLS_MAKEWINDOWS } from '../../../../modules/nf-core/bedtools/makewindows/main.nf'

workflow test_bam_ngscheckmate_bam {

    inputBed = [ [ id:'test_bed'],
                file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
                ]

    input   = [[
            [ id:'test1' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
        ],[
            [ id:'test2' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam'], checkIfExists: true)
        ]
    ]

    fasta    = [ [ id:'sarscov2'],
                file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
                ]

    BEDTOOLS_MAKEWINDOWS( Channel.of(inputBed) )
    BAM_NGSCHECKMATE(Channel.fromList(input), BEDTOOLS_MAKEWINDOWS.out.bed , Channel.of(fasta))

}

workflow test_bam_ngscheckmate_cram {

    input   = [[
            [ id:'test1' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true)
        ],[
            [ id:'test2' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_cram'], checkIfExists: true)
        ]
    ]

    inputBed  = [ [ id:'snp_bed'],
                file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_bed'], checkIfExists: true)
                ]

    fasta    = [ [ id:'homo_sapiens'],
                file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
                ]

    BEDTOOLS_MAKEWINDOWS( Channel.of(inputBed) )
    BAM_NGSCHECKMATE(Channel.fromList(input), BEDTOOLS_MAKEWINDOWS.out.bed , Channel.of(fasta))

}
