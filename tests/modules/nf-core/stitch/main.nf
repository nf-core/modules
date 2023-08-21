#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STITCH as STITCH_ONE_STEP        } from '../../../../modules/nf-core/stitch/main.nf'
include { STITCH as STITCH_GENERATE_INPUTS } from '../../../../modules/nf-core/stitch/main.nf'
include { STITCH as STITCH_IMPUTE_ONLY     } from '../../../../modules/nf-core/stitch/main.nf'

// positions and essential parameters
def posfile             = file(params.test_data['homo_sapiens']['genome']['genome_21_stitch_posfile'], checkIfExists: true)
def input_empty         = []
def rdata_empty         = []
def chromosome_name_val = "chr21"
def K_val               = 2
def nGen_val            = 1
def stitch_input        = [ [ id: "test_positions" ], posfile, input_empty, rdata_empty, chromosome_name_val, K_val, nGen_val ]

// sequencing data in cram format
def crams_val = [
    params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram' ],
    params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'],
]
def crais_val = [
    params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram_crai' ],
    params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram_crai'],
]
def reads = [ [ id:"test_reads" ], crams_val, crais_val ]

// sequencing data in bam format
def bams_val = [
    params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam' ],
    params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'],
]
def bais_val = [
    params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai' ],
    params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'],
]
def reads_bam = [ [ id:"test_reads" ], bams_val, bais_val ]

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
    cramlist = Channel.fromPath( crams_val )
    .map { it[-1] as String } // get only filename
    .collectFile( name: "cramlist.txt", newLine: true, sort: true )


    emit:
    Channel.of( reads ).combine( cramlist ).first()
}

workflow GET_READS_BAM {
    main:
    bamlist = Channel.fromPath( bams_val )
    .map { it[-1] as String } // get only filename
    .collectFile( name: "bamlist.txt", newLine: true, sort: true )


    emit:
    Channel.of( reads ).combine( bamlist ).first()
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

    Channel.of( stitch_input )
    .map {
        meta, positions, input, rdata, chromosome_name, K, nGen ->
        [ meta, positions ]
    }
    .join ( STITCH_GENERATE_INPUTS.out.input )
    .join ( STITCH_GENERATE_INPUTS.out.rdata )
    .map {
        meta, positions, input, rdata ->
        [ meta, positions, input, rdata, chromosome_name_val, K_val, nGen_val ]
    }
    .set { stitch_input_second_step }

    STITCH_IMPUTE_ONLY(
        stitch_input_second_step,
        [[id: null], [], [], []],
        [[id: null], [], []],
        seed,
    )
}

workflow test_bam {
    GET_READS_BAM()
    STITCH_ONE_STEP (
        stitch_input,
        GET_READS_BAM.out,
        reference,
        seed,
    )
}
