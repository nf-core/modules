#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STITCH } from '../../../../modules/nf-core/stitch/main.nf'

// positions and essential parameters
def posfile         = file(params.test_data['homo_sapiens']['genome']['genome_21_stitch_posfile'], checkIfExists: true)
def input           = []
def rdata           = []
def chromosome_name = "chr21"
def K               = 2
def nGen            = 1
def stitch_input    = [ [ id: "test_positions" ], input, rdata, chromosome_name, K, nGen ]

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


workflow test_stitch {

    cramlist = Channel.fromPath( crams )
    .map { it[-1] as String } // get only filename
    .collectFile( name: "cramlist.txt", newLine: true, sort: true )

    STITCH (
        stitch_input,
        reads.combine ( cramlist ),
        reference,
        seed
    )
}
