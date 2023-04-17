#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKGTF } from '../../../../../modules/nf-core/cellranger/mkgtf/main.nf'
include { CELLRANGER_MKREF } from '../../../../../modules/nf-core/cellranger/mkref/main.nf'
include { CELLRANGER_MULTI } from '../../../../../modules/nf-core/cellranger/multi/main.nf'

// cellranger multi reads everything from a "CSV"-like config
// will manually write this config and stage the corresponding files
csv_config = "cellranger_multi_config.csv"


// stage B cell FASTQ test data
bcell_fastq_1 = file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_vdj/subsampled_sc5p_v2_hs_B_1k_b_fastqs/subsampled_sc5p_v2_hs_B_1k_b_S1_L001_R1_001.fastq.gz", checkIfExists: true)
bcell_fastq_2 = file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_vdj/subsampled_sc5p_v2_hs_B_1k_b_fastqs/subsampled_sc5p_v2_hs_B_1k_b_S1_L001_R2_001.fastq.gz", checkIfExists: true)

bcell_localdir = "${workDir}/bcell_fastqs"
bcell_fastq_localdir = file( bcell_localdir )
bcell_fastq_localdir.mkdirs()

bcell_fastq_1_localpath = "${bcell_localdir}/bcells_S1_L001_R1_001.fastq.gz"
bcell_fastq_2_localpath = "${bcell_localdir}/bcells_S1_L001_R2_001.fastq.gz"

bcell_fastq_1.copyTo( bcell_fastq_1_localpath )
bcell_fastq_2.copyTo( bcell_fastq_2_localpath )


// stage 5' gene expression FASTQ test data
fivepgex_fastq_1 = file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_multi/subsampled_sc5p_v2_hs_B_1k_5gex_fastqs/subsampled_sc5p_v2_hs_B_1k_5gex_S1_L001_R1_001.fastq.gz", checkIfExists: true)
fivepgex_fastq_2 = file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_multi/subsampled_sc5p_v2_hs_B_1k_5gex_fastqs/subsampled_sc5p_v2_hs_B_1k_5gex_S1_L001_R2_001.fastq.gz", checkIfExists: true)

fivepgex_localdir = "${workDir}/5gex_fastqs"
fivepgex_fastq_localdir = file( fivepgex_localdir )
fivepgex_fastq_localdir.mkdirs()

fivepgex_fastq_1_localpath = "${fivepgex_localdir}/5gex_S1_L001_R1_001.fastq.gz"
fivepgex_fastq_2_localpath = "${fivepgex_localdir}/5gex_S1_L001_R2_001.fastq.gz"

fivepgex_fastq_1.copyTo( fivepgex_fastq_1_localpath )
fivepgex_fastq_2.copyTo( fivepgex_fastq_2_localpath )

// stage GEX reference
// will build this as done in cellranger count test
gex_fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
gex_gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
gex_reference_name = "homo_sapiens_chr22_reference"

// stage VDJ reference
vdj_reference_json = file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0/reference.json", checkIfExists: true)
vdj_reference_fasta = file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0/fasta/regions.fa", checkIfExists: true)
vdj_reference_suppfasta = file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/10xgenomics/cellranger_vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0/fasta/supp_regions.fa", checkIfExists: true)

vdj_reference_localdir = "${workDir}/vdj_reference"
vdj_reference = file( vdj_reference_localdir )
vdj_reference.mkdirs()
vdj_reference_json.copyTo("${vdj_reference_localdir}/reference.json")
vdj_reference_fasta.copyTo("${vdj_reference_localdir}/fasta/regions.fa")
vdj_reference_suppfasta.copyTo("${vdj_reference_localdir}/fasta/supp_regions.fa") 

workflow test_cellranger_multi {
    
    CELLRANGER_MKGTF ( gex_gtf )

    CELLRANGER_MKREF (
        gex_fasta,
        CELLRANGER_MKGTF.out.gtf,
        gex_reference_name
    )

    // subworkflow to create cellranger multi config file
    // also mixes reference, fastq channels for multi input
    PREPARE_MULTI_RUN(
        csv_config,
        CELLRANGER_MKREF.out.reference,
        vdj_reference,
        fivepgex_fastq_localdir,
        bcell_fastq_localdir
    )

    CELLRANGER_MULTI (
        PREPARE_MULTI_RUN.config,
        PREPARE_MULTI_RUN.references,
        PREPARE_MULTI_RUN.fastqs
    )

}


process WRITE_CONFIG {

    input:
    val config_name
    path gex_ref
    path vdj_ref
    path gex_fastqs
    path bcell_fastqs

    output:
    path("cellranger_multi_config.csv"), emit: config

    script:
    def gex_ref_dir = gex_ref.name
    def vdj_dir = vdj_ref.name
    def gex_dir = gex_fastqs.name
    def bcell_dir = bcell_fastqs.name
    def workdir = task.workDir

    """
    config="cellranger_multi_config.csv"
    touch \$config 
	echo '[gene-expression]' > \$config 
    echo 'reference,$gex_ref_dir' >> \$config
	echo 'expect-cells,1000' >> \$config
	echo 'chemistry,SC5P-PE' >> \$config
	echo '\n[vdj]' >> \$config
	echo 'reference,$vdj_dir' >> \$config
	echo '\n[libraries]' >> \$config
	echo 'fastq_id,fastqs,lanes,feature_types,subsample_rate' >> \$config
	echo '5gex,$gex_dir,1,gene expression,' >> \$config
	echo 'bcells,$bcell_dir,1,vdj,' >> \$config
    """

}

workflow PREPARE_MULTI_RUN {

    take:
    config_name
    gex_ref
    vdj_ref
    gex_fastqs
    bcell_fastqs

    main:
    WRITE_CONFIG (
        config_name,
        gex_ref,
        vdj_reference,
        gex_fastqs,
        bcell_fastqs
    )
    .out
    .config
    .map { configfile -> [[ id:'test', single_end:false ], configfile ] }
    .set { multi_input }

    emit:
    config     = multi_input
    references = Channel.mix(gex_ref, vdj_ref )
    fastqs     = Channel.mix( gex_fastqs, bcell_fastqs )
}
