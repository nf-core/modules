#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMCOPY_MAPCOUNTER } from '../../../../modules/hmmcopy/mapcounter/main.nf'
include { HMMCOPY_GENERATEMAP } from '../../../../modules/hmmcopy/generatemap/main.nf'
workflow test_hmmcopy_mapcounter {

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    HMMCOPY_GENERATEMAP( fasta )

    HMMCOPY_MAPCOUNTER ( HMMCOPY_GENERATEMAP.out.bigwig )
}
