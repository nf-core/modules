#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STITCH as STITCH_ONE_STEP        } from '../../../../modules/nf-core/stitch/main.nf'
include { STITCH as STITCH_GENERATE_INPUTS } from '../../../../modules/nf-core/stitch/main.nf'
include { STITCH as STITCH_IMPUTE_ONLY     } from '../../../../modules/nf-core/stitch/main.nf'

// positions and essential parameters
def posfile         = file(params.test_data['homo_sapiens']['genome']['genome_21_stitch_posfile'], checkIfExists: true)
def input           = []
def rdata           = []
def chromosome_name = "chr21"
def K               = 2
def nGen            = 1
def stitch_input    = [ [ id: "test_positions" ], posfile, input, rdata, chromosome_name, K, nGen ]

// sequencing data
def crams = [
    params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram' ],
    params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'],
]
def crais = [
    params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram_crai' ],
    params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram_crai'],
]
def reads = [ [ id:"test_reads" ], crams, crais ]

// reference genome
def reference = [
    [ id:"test_reference" ],
    file(params.test_data['homo_sapiens']['genome']['genome_21_fasta']    , checkIfExists: true),
    file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true),
]

// for reproducibility
def seed = 1

//
// Test workflows
//

workflow GET_READS {
    main:
    cramlist = Channel.fromPath( crams )
    .map { it[-1] as String } // get only filename
    .collectFile( name: "cramlist.txt", newLine: true, sort: true )


    emit:
    Channel.of( reads ).combine( cramlist ).first()
}


workflow test_with_seed {
    GET_READS()
    STITCH_ONE_STEP (
        stitch_input,
        GET_READS.out,
        reference,
        seed,
    )
}

workflow test_no_seed {
    GET_READS()
    STITCH_ONE_STEP (
        stitch_input,
        GET_READS.out,
        reference,
        [],
    )
}

workflow test_two_stage_imputation {
    GET_READS()
    STITCH_GENERATE_INPUTS (
        stitch_input,
        GET_READS.out,
        reference,
        seed,
    )

    stitch_input
    .map {
        meta, positions, input, rdata, chromosome_name, K, nGen ->
        [ meta, positions, chromosome_name ]
    }
    .join ( STITCH_GENERATE_INPUTS.out.input )
    .join ( STITCH_GENERATE_INPUTS.out.rdata )
    .map {
        meta, positions, chromosome_name, input, rdata ->
        [ meta, positions, input, rdata, chromosome_name ]
    }
    .combine ( K    )
    .combine ( nGen )
    .set { stitch_input_second_step }

    STITCH_IMPUTE_ONLY(
        stitch_input_second_step,
        [[id: null], [], []],
        [[id: null], [], []],
        seed,
    )
}
