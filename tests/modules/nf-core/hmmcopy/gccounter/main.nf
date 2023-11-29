#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMCOPY_GCCOUNTER } from '../../../../../modules/nf-core/hmmcopy/gccounter/main.nf'

workflow test_hmmcopy_gccounter {
    fasta = Channel.of(
        tuple(
            [id: "test"],
            file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        )
    )

    HMMCOPY_GCCOUNTER (fasta)
}
