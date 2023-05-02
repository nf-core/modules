#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKGTF } from '../../../../../modules/nf-core/cellranger/mkgtf/main.nf'
include { CELLRANGER_MKREF } from '../../../../../modules/nf-core/cellranger/mkref/main.nf'
include { CELLRANGER_MULTI } from '../../../../../modules/nf-core/cellranger/multi/main.nf'

// stage B cell FASTQ test data
bcell_fastqs = [
    file(params.test_data['homo_sapiens']['illumina']['test_10x_b_1_fastq_gz'], checkIfExists: true),
    file(params.test_data['homo_sapiens']['illumina']['test_10x_b_2_fastq_gz'], checkIfExists: true)
]
def bcell_fastq_samplename = "subsampled_sc5p_v2_hs_B_1k_b"

// stage 5' gene expression FASTQ test data
gex_fastqs = [
    file(params.test_data['homo_sapiens']['illumina']['test_10x_5p_1_fastq_gz'], checkIfExists: true),
    file(params.test_data['homo_sapiens']['illumina']['test_10x_5p_2_fastq_gz'], checkIfExists: true)
]
def gex_fastq_samplename = "subsampled_sc5p_v2_hs_B_1k_5gex"

// stage GEX reference
// will build this as done in cellranger count test
gex_ref_fasta    = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
gex_ref_gtf      = file(params.test_data['homo_sapiens']['genome']['genome_gtf']  , checkIfExists: true)
def gex_ref_name = "homo_sapiens_chr22_reference"

// stage VDJ reference
vdj_json      = file(params.test_data['homo_sapiens']['illumina']['test_10x_vdj_ref_json']     , checkIfExists: true)
vdj_fasta     = file(params.test_data['homo_sapiens']['illumina']['test_10x_vdj_ref_fasta']    , checkIfExists: true)
vdj_suppfasta = file(params.test_data['homo_sapiens']['illumina']['test_10x_vdj_ref_suppfasta'], checkIfExists: true)

// awkwardly restage VDJ ref to enforce directory structure expected by cellranger
def vdj_ref_name      = "vdj_reference"
def vdj_reference_dir = "$workDir/$vdj_ref_name"
vdj_reference         = file( vdj_reference_dir )
vdj_json.copyTo("$vdj_reference_dir/reference.json")
vdj_fasta.copyTo("$vdj_reference_dir/fasta/regions.fa")
vdj_suppfasta.copyTo("$vdj_reference_dir/fasta/supp_regions.fa")

// make an empty dummy file
empty_file = file("$workDir/EMPTY")
empty_file.append("")


// create empty channels to fill unused cellranger multi arguments
// fastqs need a [ meta, ref ] structure
// references just need a path
ch_ab_fastqs             = Channel.fromPath( empty_file ).map { file -> [ [ id:"EMPTY", options:[] ], file ] }
ch_beam_fastqs           = Channel.fromPath( empty_file ).map { file -> [ [ id:"EMPTY", options:[] ], file ] }
ch_cmo_fastqs            = Channel.fromPath( empty_file ).map { file -> [ [ id:"EMPTY", options:[] ], file ] }
ch_crispr_fastqs         = Channel.fromPath( empty_file ).map { file -> [ [ id:"EMPTY", options:[] ], file ] }
ch_gex_frna_probeset     = Channel.fromPath( empty_file )
ch_gex_targetpanel       = Channel.fromPath( empty_file )
ch_vdj_primer_index      = Channel.fromPath( empty_file )
ch_fb_reference          = Channel.fromPath( empty_file )
ch_beam_panel            = Channel.fromPath( empty_file )
ch_cmo_reference         = Channel.fromPath( empty_file )
ch_cmo_barcodes          = Channel.fromPath( empty_file )
ch_cmo_sample_assignment = Channel.fromPath( empty_file )
ch_frna_sampleinfo       = Channel.fromPath( empty_file )

def test_meta = [ id:'test', single_end:false ]


workflow test_cellranger_multi {

    CELLRANGER_MKGTF ( gex_ref_gtf )

    CELLRANGER_MKREF ( gex_ref_fasta, CELLRANGER_MKGTF.out.gtf, gex_ref_name)

    // collect references and fastq files for staging
    ch_gex_reference = CELLRANGER_MKREF.out.reference
    ch_gex_fastqs    = Channel.of( gex_fastqs )
        .collect()
        .map { reads -> [ [ id:gex_fastq_samplename, options:[ "expect-cells":"1000", chemistry:"SC5P-PE" ] ], reads ] }
    ch_vdj_reference = Channel.of( vdj_reference )
    ch_bcell_fastqs  = Channel.of( bcell_fastqs )
        .collect()
        .map { reads -> [ [ id:bcell_fastq_samplename, options:[] ], reads ] }


    CELLRANGER_MULTI (
        test_meta,
        ch_gex_fastqs,
        ch_bcell_fastqs,
        ch_ab_fastqs,
        ch_beam_fastqs,
        ch_cmo_fastqs,
        ch_crispr_fastqs,
        CELLRANGER_MKREF.out.reference,
        ch_gex_frna_probeset,
        ch_gex_targetpanel,
        ch_vdj_reference,
        ch_vdj_primer_index,
        ch_fb_reference,
        ch_beam_panel,
        ch_cmo_reference,
        ch_cmo_barcodes,
        ch_cmo_sample_assignment,
        ch_frna_sampleinfo
    )

}
