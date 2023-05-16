#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKGTF } from '../../../../../modules/nf-core/cellranger/mkgtf/main.nf'
include { CELLRANGER_MKREF } from '../../../../../modules/nf-core/cellranger/mkref/main.nf'
include { CELLRANGER_MULTI } from '../../../../../modules/nf-core/cellranger/multi/main.nf'

/***************************/
/*** stage 10k PBMC data ***/
/***************************/

// stage B cell FASTQ test data
bcell_fastqs_10k_pbmc = [
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_b_fastq_1_gz'], checkIfExists: true),
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_b_fastq_2_gz'], checkIfExists: true)
]
def bcell_fastq_samplename_10k_pbmc = "subsampled_sc5p_v2_hs_PBMC_10k_b"

// stage 5' gene expression FASTQ test data
fivepgex_fastqs_10k_pbmc = [
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_5gex_fastq_1_gz'], checkIfExists: true),
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_5gex_fastq_2_gz'], checkIfExists: true)
]
def fivepgex_fastq_samplename_10k_pbmc = "subsampled_sc5p_v2_hs_PBMC_10k_5gex"

// stage 5' feature barcode (antibody capture) FASTQ test data
ab_fastqs_10k_pbmc = [
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_5fb_fastq_1_gz'], checkIfExists: true),
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_5fb_fastq_2_gz'], checkIfExists: true)
]
def fivepab_fastq_samplename_10k_pbmc = "subsampled_sc5p_v2_hs_PBMC_10k_5fb"

// stage feature barcode reference for antibody capture
fb_reference_10k_pbmc = file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_feature_ref_csv'], checkIfExists: true)

/*******************************/
/*** end stage 10k PBMC data ***/
/*******************************/

/**********************************/
/*** stage 10k PBMC w/ CMO data ***/
/**********************************/

//test_10x_10k_pbmc_cmo_gex2_fastq_1_gz
//test_10x_10k_pbmc_cmo_gex2_fastq_2_gz

// stage 3' CMO FASTQ test data
cmo_fastqs_10k_pbmc_cmo = [
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_cmo_cmo_fastq_1_gz'], checkIfExists: true),
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_cmo_cmo_fastq_2_gz'], checkIfExists: true)
]
def cmo_fastq_samplename_10k_pbmc_cmo = "subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture"

// stage 3' gene expression FASTQ test data
threepgex_fastqs_10k_pbmc_cmo = [
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_cmo_gex1_fastq_1_gz'], checkIfExists: true),
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_cmo_gex1_fastq_2_gz'], checkIfExists: true)
]
def threepgex_fastq_samplename_10k_pbmc_cmo = "subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex"

// stage feature barcode reference for antibody capture
cmo_reference_10k_pbmc_cmo = file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_cmo_feature_ref_csv'], checkIfExists: true)

// CMO needs a barcode file
cmo_barcodes_csv = file("cmo_barcodes.csv")
cmo_barcodes_csv.text = ""
cmo_barcodes_csv.append("sample_id,cmo_ids,description\nPBMCs_human_1,CMO301,PBMCs_human_1\nPBMCs_human_2,CMO302,PBMCs_human_2")


/**************************************/
/*** end stage 10k PBMC w/ CMO data ***/
/**************************************/



/*********************************/
/*** stage 5k CMV+ T-cell data ***/
/*********************************/

// stage antibody capture data 
ab_fastqs_5k_cmvpos_tcells = [
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_ab_fastq_1_gz'], checkIfExists: true),
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_ab_fastq_2_gz'], checkIfExists: true)
]
def ab_fastq_samplename_5k_cmvpos_tcells = "subsampled_5k_human_antiCMV_T_TBNK_connect_AB"

// stage GEX data 
gex_fastqs_5k_cmvpos_tcells = [
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_gex1_fastq_1_gz'], checkIfExists: true),
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_gex1_fastq_2_gz'], checkIfExists: true)
]
def gex_fastq_samplename_5k_cmvpos_tcells = "subsampled_5k_human_antiCMV_T_TBNK_connect_GEX_1"

// stage VDJ data 
vdj_fastqs_5k_cmvpos_tcells = [ 
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_vdj_fastq_1_gz'], checkIfExists: true),
    file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_vdj_fastq_2_gz'], checkIfExists: true)
]
def vdj_fastq_samplename_5k_cmvpos_tcells = "subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ"

// stage feature barcode reference for antibody capture
fb_reference_5k_cmvpos_tcells = file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_feature_ref_csv'], checkIfExists: true)

/*************************************/
/*** end stage 5k CMV+ T-cell data ***/
/*************************************/



/***************************/
/*** stage GEX reference ***/
/***************************/

// will build this as done in cellranger count test
gex_ref_fasta    = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
gex_ref_gtf      = file(params.test_data['homo_sapiens']['genome']['genome_gtf']  , checkIfExists: true)
def gex_ref_name = "homo_sapiens_chr22_reference"

/*******************************/
/*** end stage GEX reference ***/
/*******************************/


/***************************/
/*** stage VDJ reference ***/
/***************************/

vdj_json      = file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_vdj_ref_json']     , checkIfExists: true)
vdj_fasta     = file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_vdj_ref_fasta']    , checkIfExists: true)
vdj_suppfasta = file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_vdj_ref_suppfasta'], checkIfExists: true)

// awkwardly restage VDJ ref to enforce directory structure expected by cellranger
def vdj_ref_name      = "vdj_reference"
def vdj_reference_dir = "$workDir/$vdj_ref_name"
vdj_reference         = file( vdj_reference_dir )
vdj_json.copyTo("$vdj_reference_dir/reference.json")
vdj_fasta.copyTo("$vdj_reference_dir/fasta/regions.fa")
vdj_suppfasta.copyTo("$vdj_reference_dir/fasta/supp_regions.fa")

/*******************************/
/*** end stage VDJ reference ***/
/*******************************/


// make an empty dummy file
empty_file = file("$workDir/EMPTY")
empty_file.append("")


// create empty channels to fill unused cellranger multi arguments
// fastqs need a [ meta, ref ] structure
// references just need a path
ch_gex_fastqs            = Channel.fromPath( empty_file ).map { file -> [ [ id:"EMPTY", options:[] ], file ] }
ch_vdj_fastqs            = Channel.fromPath( empty_file ).map { file -> [ [ id:"EMPTY", options:[] ], file ] }
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


workflow test_cellranger_multi_10k_pbmc {

    def test_meta_10k_pbmc = [ id:'test_10k_pbmc', single_end:false ]

    // collect references and fastq files for staging
    ch_gex_fastqs_10k_pbmc = Channel.of( fivepgex_fastqs_10k_pbmc )
        .collect()
        .map { reads -> [ [ id:fivepgex_fastq_samplename_10k_pbmc , options:[ "expect-cells":"1000", chemistry:"SC5P-PE", "no-bam":true, "no-secondary":true ] ], reads ] }

    ch_vdj_reference = Channel.of( vdj_reference )

    ch_bcell_fastqs_10k_pbmc  = Channel.of( bcell_fastqs_10k_pbmc )
        .collect()
        .map { reads -> [ [ id:bcell_fastq_samplename_10k_pbmc , options:[] ], reads ] }

    ch_ab_fastqs_10k_pbmc  = Channel.of( ab_fastqs_10k_pbmc )
        .collect()
        .map { reads -> [ [ id:fivepab_fastq_samplename_10k_pbmc  , options:[] ], reads ] }

    ch_ab_reference_10k_pbmc = Channel.of( fb_reference_10k_pbmc )

    CELLRANGER_MKGTF ( gex_ref_gtf )

    CELLRANGER_MKREF ( gex_ref_fasta, CELLRANGER_MKGTF.out.gtf, gex_ref_name)

    CELLRANGER_MULTI (
        test_meta_10k_pbmc,
        ch_gex_fastqs_10k_pbmc,
        ch_bcell_fastqs_10k_pbmc,
        ch_ab_fastqs_10k_pbmc,
        ch_beam_fastqs,
        ch_cmo_fastqs,
        ch_crispr_fastqs,
        CELLRANGER_MKREF.out.reference,
        ch_gex_frna_probeset,
        ch_gex_targetpanel,
        ch_vdj_reference,
        ch_vdj_primer_index,
        ch_ab_reference_10k_pbmc,
        ch_beam_panel,
        ch_cmo_reference,
        ch_cmo_barcodes,
        ch_cmo_sample_assignment,
        ch_frna_sampleinfo
    )

}

workflow test_cellranger_multi_10k_pbmc_cmo {

    def test_meta_10k_pbmc_cmo = [ id:'test_10k_pbmc_cmo', single_end:false ]

    // collect references and fastq files for staging
    ch_gex_fastqs_10k_pbmc_cmo = Channel.of( threepgex_fastqs_10k_pbmc_cmo )
        .collect()
        .map { reads -> [ [ id:threepgex_fastq_samplename_10k_pbmc_cmo , options:[ "expect-cells":"1000", chemistry:"SC3Pv3", "no-bam":true, "no-secondary":true ] ], reads ] }

    ch_vdj_reference = Channel.of( vdj_reference )

    ch_cmo_fastqs_10k_pbmc_cmo  = Channel.of( cmo_fastqs_10k_pbmc_cmo )
        .collect()
        .map { reads -> [ [ id:cmo_fastq_samplename_10k_pbmc_cmo, options:[] ], reads ] }

    ch_cmo_reference_10k_pbmc_cmo = Channel.of( cmo_reference_10k_pbmc_cmo )

    // empty vdj reference
    ch_vdj_ref_empty = Channel.fromPath( empty_file ) 

    // CMO analysis needs barcodes
    ch_cmo_barcodes_10k_pbmc_cmo = Channel.fromPath( cmo_barcodes_csv )

    CELLRANGER_MKGTF ( gex_ref_gtf )

    CELLRANGER_MKREF ( gex_ref_fasta, CELLRANGER_MKGTF.out.gtf, gex_ref_name)

    CELLRANGER_MULTI (
        test_meta_10k_pbmc_cmo,
        ch_gex_fastqs_10k_pbmc_cmo,
        ch_vdj_fastqs,
        ch_ab_fastqs,
        ch_beam_fastqs,
        ch_cmo_fastqs_10k_pbmc_cmo,
        ch_crispr_fastqs,
        CELLRANGER_MKREF.out.reference,
        ch_gex_frna_probeset,
        ch_gex_targetpanel,
        ch_vdj_ref_empty,
        ch_vdj_primer_index,
        ch_fb_reference,
        ch_beam_panel,
        ch_cmo_reference_10k_pbmc_cmo,
        ch_cmo_barcodes_10k_pbmc_cmo,
        ch_cmo_sample_assignment,
        ch_frna_sampleinfo
    )
}

workflow test_cellranger_multi_5k_cmvpos_tcells {

    def test_meta_5k_cmvpos_tcells = [ id:'test_5k_cmvpos_tcells', single_end:false ]

    // collect references and fastq files for staging
    ch_gex_fastqs_5k_cmvpos_tcells = Channel.of( gex_fastqs_5k_cmvpos_tcells )
        .collect()
        .map { reads -> [ [ id:gex_fastq_samplename_5k_cmvpos_tcells , options:[ "expect-cells":"1000", chemistry:"SC5P-R2", "no-bam":true, "no-secondary":true ] ], reads ] }

    ch_vdj_reference = Channel.of( vdj_reference )

    ch_vdj_fastqs_5k_cmvpos_tcells  = Channel.of( vdj_fastqs_5k_cmvpos_tcells )
        .collect()
        .map { reads -> [ [ id:vdj_fastq_samplename_5k_cmvpos_tcells, options:[] ], reads ] }

    ch_fb_reference_5k_cmvpos_tcells = Channel.of( fb_reference_5k_cmvpos_tcells )

    ch_ab_fastqs_5k_cmvpos_tcells  = Channel.of( ab_fastqs_5k_cmvpos_tcells )
        .collect()
        .map { reads -> [ [ id:ab_fastq_samplename_5k_cmvpos_tcells, options:[] ], reads ] }

    CELLRANGER_MKGTF ( gex_ref_gtf )

    CELLRANGER_MKREF ( gex_ref_fasta, CELLRANGER_MKGTF.out.gtf, gex_ref_name)

    CELLRANGER_MULTI (
        test_meta_5k_cmvpos_tcells,
        ch_gex_fastqs_5k_cmvpos_tcells,
        ch_vdj_fastqs_5k_cmvpos_tcells,
        ch_ab_fastqs_5k_cmvpos_tcells,
        ch_beam_fastqs,
        ch_cmo_fastqs,
        ch_crispr_fastqs,
        CELLRANGER_MKREF.out.reference,
        ch_gex_frna_probeset,
        ch_gex_targetpanel,
        ch_vdj_reference,
        ch_vdj_primer_index,
        ch_fb_reference_5k_cmvpos_tcells,
        ch_beam_panel,
        ch_cmo_reference,
        ch_cmo_barcodes,
        ch_cmo_sample_assignment,
        ch_frna_sampleinfo
    )
}
