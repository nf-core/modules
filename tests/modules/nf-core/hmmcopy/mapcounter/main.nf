#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMCOPY_MAPCOUNTER } from '../../../../../modules/nf-core/hmmcopy/mapcounter/main.nf'
include { HMMCOPY_GENERATEMAP } from '../../../../../modules/nf-core/hmmcopy/generatemap/main.nf'
workflow test_hmmcopy_mapcounter {

    fasta = Channel.of(
        tuple(
            [id: "test"],
            file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        )
    )

    HMMCOPY_GENERATEMAP( fasta )

    HMMCOPY_MAPCOUNTER ( HMMCOPY_GENERATEMAP.out.bigwig )
}
