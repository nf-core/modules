#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HLALA_PREPAREGRAPH } from '../../../../../modules/nf-core/hlala/preparegraph/main.nf'

workflow test_hlala_preparegraph {
    HLALA_PREPAREGRAPH ( Channel.fromPath("/home-link/iivow01/git/modules/hla_graph/PRG_MHC_GRCh38_withIMGT") )
}
