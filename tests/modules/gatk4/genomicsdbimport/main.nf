#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR           } from '../../../../modules/untar/main.nf'           addParams( options: [:] )
include { GATK4_GENOMICSDBIMPORT } from '../../../../modules/gatk4/genomicsdbimport/main.nf' addParams( options: [:] )

workflow test_gatk4_genomicsdbimport_create_genomicsdb {

    input = [ [ id:'test'], // meta map
              file(params.test_data['homo_sapiens']['genome']['syntheticvcf_short_vcf_gz'], checkIfExists: true) ,
              file(params.test_data['homo_sapiens']['genome']['syntheticvcf_short_vcf_gz_tbi'], checkIfExists: true) ,
              [] ,
              '22' ,
              [] ]

    run_intlist = false
    run_updatewspace = false
    input_map = false

    GATK4_GENOMICSDBIMPORT ( input, run_intlist, run_updatewspace, input_map )
}

workflow test_gatk4_genomicsdbimport_get_intervalslist {
    db    = file(params.test_data['homo_sapiens']['illumina']['test_genomicsdb_tar_gz'], checkIfExists: true)

    UNTAR ( db )

    def input = Channel.of([ [ id:'test'], // meta map
              [] ,
              [] ,
              [] ,
              [] ])
              .combine(UNTAR.out.untar)

    run_intlist = true
    run_updatewspace = false
    input_map = false

    GATK4_GENOMICSDBIMPORT ( input, run_intlist, run_updatewspace, input_map )
}

workflow test_gatk4_genomicsdbimport_update_genomicsdb {
    db    = file(params.test_data['homo_sapiens']['illumina']['test_genomicsdb_tar_gz'], checkIfExists: true)

    UNTAR ( db )

    def input = Channel.of([ [ id:'test'], // meta map
              file( params.test_data['homo_sapiens']['illumina']['test2_genome_vcf_gz'] , checkIfExists: true) ,
              file( params.test_data['homo_sapiens']['illumina']['test2_genome_vcf_gz_tbi'] , checkIfExists: true) ,
              [] ,
              [] ])
              .combine(UNTAR.out.untar)

    run_intlist = false
    run_updatewspace = true
    input_map = false

    GATK4_GENOMICSDBIMPORT ( input, run_intlist, run_updatewspace, input_map )

}
