#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENOMAD_DOWNLOAD } from '../../../../../modules/nf-core/genomad/download/main.nf'
include { GENOMAD_ENDTOEND } from '../../../../../modules/nf-core/genomad/endtoend/main.nf'

workflow test_genomad_endtoend {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('https://github.com/nf-core/test-datasets/raw/mag/assemblies/SPAdes-test_minigut_contigs.fasta.gz', checkIfExists: true)
    ]

    GENOMAD_DOWNLOAD ( )

    GENOMAD_ENDTOEND ( input, GENOMAD_DOWNLOAD.out.genomad_db )
}
