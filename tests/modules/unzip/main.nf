#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNZIP } from '../../../modules/unzip/main.nf' addParams( options: [:] )

workflow test_unzip {

    archive = [ file(params.test_data['sarscov2']['genome']['genome_fasta_zip'], checkIfExists: true),
                file('/home/jfellows/Downloads/ss_imgs/participants-affiliation.zip', checkIfExists: true)
                ]

    UNZIP ( archive )
}
