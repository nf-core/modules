#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOBSUITE_RECON } from '../../../../modules/mobsuite/recon/main.nf'

workflow test_mobsuite_recon {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    MOBSUITE_RECON ( input )
}
