#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_GENOMICSDBIMPORT } from '../../../../modules/gatk4/genomicsdbimport/main.nf' addParams( options: [:] )

workflow test_gatk4_genomicsdbimport_create_genomicsdb {

    input = [ [ id:'test'], // meta map
              file( params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'] , checkIfExists: true) ,
              file( params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'] , checkIfExists: true) ,
              [] ,
              file( "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.interval_list" , checkIfExists: true) ,
              [] ]

    run_intlist = false
    run_updatewspace = false
    input_map = false

    GATK4_GENOMICSDBIMPORT ( input , run_intlist , run_updatewspace , input_map )
}
