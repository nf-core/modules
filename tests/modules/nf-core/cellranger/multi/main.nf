#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKGTF } from '../../../../../modules/nf-core/cellranger/mkgtf/main.nf'
include { CELLRANGER_MKREF } from '../../../../../modules/nf-core/cellranger/mkref/main.nf'
include { CELLRANGER_MULTI } from '../../../../../modules/nf-core/cellranger/multi/main.nf'

// cellranger multi reads everything from a "CSV"-like config
// will manually write this config and stage the corresponding files
csv_config = "cellranger_multi_config.csv"


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
gex_ref_fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
gex_ref_gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
def gex_ref_name = "homo_sapiens_chr22_reference"

// stage VDJ reference
// must stage awkwardly to preserve directory structure expected by cellranger
vdj_json = file(params.test_data['homo_sapiens']['illumina']['test_10x_vdj_ref_json'], checkIfExists: true)
vdj_fasta = file(params.test_data['homo_sapiens']['illumina']['test_10x_vdj_ref_fasta'], checkIfExists: true)
vdj_suppfasta = file(params.test_data['homo_sapiens']['illumina']['test_10x_vdj_ref_suppfasta'], checkIfExists: true)

def vdj_ref_name = "vdj_reference"
def vdj_reference_dir = "$workDir/$vdj_ref_name"
vdj_reference = file( vdj_reference_dir )
vdj_json.copyTo("$vdj_reference_dir/reference.json")
vdj_fasta.copyTo("$vdj_reference_dir/fasta/regions.fasta")
vdj_suppfasta.copyTo("$vdj_reference_dir/fasta/supp_regions.fasta")


workflow test_cellranger_multi {
    
    CELLRANGER_MKGTF ( gex_ref_gtf )


    CELLRANGER_MKREF (
        gex_ref_fasta,
        CELLRANGER_MKGTF.out.gtf,
        gex_ref_name
    )

    // collect references and fastq files for staging
    ch_references = CELLRANGER_MKREF.out.reference
        .mix ( Channel.of ( vdj_reference ) )
        .collect()
    ch_fastqs = Channel.of( gex_fastqs )
        .mix( Channel.of( bcell_fastqs ) )
        .collect()


    WRITE_CONFIG ( csv_config )
    multi_input = WRITE_CONFIG
        .out
        .config
        .map { configfile -> [[ id:'test', single_end:false ], configfile ] }


    CELLRANGER_MULTI (
        multi_input,
        ch_references,
        ch_fastqs
    )

}


process WRITE_CONFIG {

    input:
    val config_name

    output:
    path "cellranger_multi_config.csv", emit: config

    script:
    """
    config="cellranger_multi_config.csv"
    touch \$config 
	echo "[gene-expression]" > \$config 
    echo "reference,references/$gex_ref_name" >> \$config
	echo "expect-cells,1000" >> \$config
	echo "chemistry,SC5P-PE" >> \$config
	echo "\n[vdj]" >> \$config
	echo "reference,references" >> \$config
	echo "\n[libraries]" >> \$config
	echo "fastq_id,fastqs,lanes,feature_types,subsample_rate" >> \$config
	echo "$gex_fastq_samplename,fastqs,1,gene expression," >> \$config
	echo "$bcell_fastq_samplename,fastqs,1,vdj," >> \$config
    """

}
