#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { UNTAR           } from "$moduleDir/modules/nf-core/untar/main.nf"
include { GATK4_CREATESOMATICPANELOFNORMALS } from "$moduleDir/modules/nf-core/gatk4/createsomaticpanelofnormals/main.nf"

workflow test_gatk4_createsomaticpanelofnormals {
    db    = [[], file(params.test_data['homo_sapiens']['illumina']['test_genomicsdb_tar_gz'], checkIfExists: true) ]

    UNTAR ( db )

    input = Channel.of([ id:'test'])
              .combine(UNTAR.out.untar.map{ it[1] })

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fastaidx = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_CREATESOMATICPANELOFNORMALS ( input, fasta, fastaidx, dict )
}
