#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_GENOMICSDBIMPORT } from '../../../../modules/gatk4/genomicsdbimport/main.nf' addParams( options: [:] )

workflow test_gatk4_genomicsdbimport_create_genomicsdb {

    input = [ [ id:'test'], // meta map
              file( params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'] , checkIfExists: true) ,
              file( params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'] , checkIfExists: true) ,
              [] ,
              file( params.test_data['homo_sapiens']['genome']['genome_interval_list'] , checkIfExists: true) ,
              [] ]

    run_intlist = false
    run_updatewspace = false
    input_map = false

    GATK4_GENOMICSDBIMPORT ( input, run_intlist, run_updatewspace, input_map )
}

workflow test_gatk4_genomicsdbimport_get_intervalslist {
    maindir = file('test_genomicsdb')
    subdir1 = file('test_genomicsdb/chr22$1$40001')
    subdir2 = file('test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448')
    subdir3 = file('test_genomicsdb/chr22$1$40001/genomicsdb_meta_dir')
    subdir2.mkdirs()
    subdir3.mkdirs()

    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/__tiledb_workspace.tdb' , checkIfExists: true).copyTo(maindir)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/vcfheader.vcf' , checkIfExists: true).copyTo(maindir)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/vidmap.json' , checkIfExists: true).copyTo(maindir)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/callset.json' , checkIfExists: true).copyTo(maindir)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/.__consolidation_lock' , checkIfExists: true).copyTo(subdir1)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__array_schema.tdb' , checkIfExists: true).copyTo(subdir1)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/genomicsdb_meta_dir/genomicsdb_column_bounds.json' , checkIfExists: true).copyTo(subdir3)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/genomicsdb_meta_dir/genomicsdb_meta_2b25a6c2-cb94-4a4a-9005-acb7c595d322.json' , checkIfExists: true).copyTo(subdir3)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/AD.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/AD_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ALT.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ALT_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/BaseQRankSum.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/DB.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/DP.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/DP_FORMAT.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/END.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ExcessHet.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/FILTER.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/FILTER_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/GQ.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/GT.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/GT_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ID.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ID_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/InbreedingCoeff.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MIN_DP.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MLEAC.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MLEAC_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MLEAF.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MLEAF_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MQRankSum.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PGT.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PGT_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PID.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PID_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PL.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PL_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PS.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/QUAL.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/RAW_MQandDP.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/REF.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/REF_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ReadPosRankSum.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/SB.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/__book_keeping.tdb.gz' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/__coords.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/__tiledb_fragment.tdb' , checkIfExists: true).copyTo(subdir2)

    input = [ [ id:'test'], // meta map
              [] ,
              [] ,
              file( maindir , checkIfExists: true) ,
              [] ,
              [] ]

    run_intlist = true
    run_updatewspace = false
    input_map = false

    GATK4_GENOMICSDBIMPORT ( input, run_intlist, run_updatewspace, input_map )
}

workflow test_gatk4_genomicsdbimport_update_genomicsdb {
    maindir = file('test_genomicsdb')
    subdir1 = file('test_genomicsdb/chr22$1$40001')
    subdir2 = file('test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448')
    subdir3 = file('test_genomicsdb/chr22$1$40001/genomicsdb_meta_dir')
    subdir2.mkdirs()
    subdir3.mkdirs()

    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/__tiledb_workspace.tdb' , checkIfExists: true).copyTo(maindir)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/vcfheader.vcf' , checkIfExists: true).copyTo(maindir)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/vidmap.json' , checkIfExists: true).copyTo(maindir)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/callset.json' , checkIfExists: true).copyTo(maindir)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/.__consolidation_lock' , checkIfExists: true).copyTo(subdir1)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__array_schema.tdb' , checkIfExists: true).copyTo(subdir1)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/genomicsdb_meta_dir/genomicsdb_column_bounds.json' , checkIfExists: true).copyTo(subdir3)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/genomicsdb_meta_dir/genomicsdb_meta_2b25a6c2-cb94-4a4a-9005-acb7c595d322.json' , checkIfExists: true).copyTo(subdir3)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/AD.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/AD_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ALT.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ALT_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/BaseQRankSum.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/DB.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/DP.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/DP_FORMAT.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/END.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ExcessHet.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/FILTER.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/FILTER_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/GQ.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/GT.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/GT_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ID.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ID_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/InbreedingCoeff.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MIN_DP.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MLEAC.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MLEAC_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MLEAF.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MLEAF_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/MQRankSum.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PGT.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PGT_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PID.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PID_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PL.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PL_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/PS.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/QUAL.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/RAW_MQandDP.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/REF.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/REF_var.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/ReadPosRankSum.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/SB.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/__book_keeping.tdb.gz' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/__coords.tdb' , checkIfExists: true).copyTo(subdir2)
    file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_genomicsdb/chr22$1$40001/__3cf81648-433d-4464-be08-23d082445c9b139814474716928_1630588248448/__tiledb_fragment.tdb' , checkIfExists: true).copyTo(subdir2)

    input = [ [ id:'test'], // meta map
              file( params.test_data['homo_sapiens']['illumina']['test2_genome_vcf_gz'] , checkIfExists: true) ,
              file( params.test_data['homo_sapiens']['illumina']['test2_genome_vcf_gz_tbi'] , checkIfExists: true) ,
              file( maindir , checkIfExists: true) ,
              [] ,
              [] ]
    run_intlist = false
    run_updatewspace = true
    input_map = false

    GATK4_GENOMICSDBIMPORT ( input, run_intlist, run_updatewspace, input_map )

}
