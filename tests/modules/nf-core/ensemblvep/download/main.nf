#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENSEMBLVEP_DOWNLOAD } from '../../../../../modules/nf-core/ensemblvep/download/main.nf'

workflow test_ensemblvep_download {
    input = [[id:"test"], "WBcel235", "caenorhabditis_elegans", "108"]

    ENSEMBLVEP_DOWNLOAD(input)
}
