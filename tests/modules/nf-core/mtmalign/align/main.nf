#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MTMALIGN_ALIGN } from '../../../../../modules/nf-core/mtmalign/align/main.nf'

workflow test_mtmalign_align {

    pdbs = [
        [ id:'test' ],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/structures/seatoxin-ref.tar.gz", checkIfExists: true)
    ]


    ch_pdbs = UNTAR ( structures ).untar
    ch_pdbs = ch_pdbs.map { meta,dir -> [file(dir).listFiles().collect()]}


    MTMALIGN_ALIGN ([] , ch_structures)
}
