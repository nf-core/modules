#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../../../modules/gatk4/createsomaticpanelofnormals/main.nf' addParams( options: [suffix:'.pon'] )

workflow test_gatk4_createsomaticpanelofnormals {

    input = [ [ id:'test' ], // meta map
              file( '/home/AD/gmackenz/test_data/test_genomicsdb' , checkIfExists: true)]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fastaidx = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_CREATESOMATICPANELOFNORMALS ( input, fasta, fastaidx, dict )
}
