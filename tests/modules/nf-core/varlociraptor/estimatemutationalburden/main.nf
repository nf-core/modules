#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_ESTIMATEMUTATIONALBURDEN } from '../../../../../modules/nf-core/varlociraptor/estimatemutationalburden/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_VCF       } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'

workflow test_varlociraptor_estimatemutationalburden {

    input = [
        [ id:'test' ], // meta map
        file("./HCC1395T_vs_HCC1395N.strelka.somatic_indels_VEP.ann.vcf.gz", checkIfExists: true),
    ]

    VARLOCIRAPTOR_ESTIMATEMUTATIONALBURDEN ( input , "hist" )
}
