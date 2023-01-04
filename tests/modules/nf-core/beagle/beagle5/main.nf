#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEAGLE_BEAGLE5 } from '../../../../../modules/nf-core/beagle/beagle5/main.nf'

workflow test_beagle_beagle5 {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/ashotmarg/modules/raw/master/ashTestData/target.22Jul22.46e.vcf.gz", checkIfExists: true)
    ]
    []
    []
    []
    []

    BEAGLE_BEAGLE5 ( input,[],[],[],[] )
}

workflow test_beagle_beagle5_ref {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/ashotmarg/modules/raw/master/ashTestData/target.22Jul22.46e.vcf.gz", checkIfExists: true)
    ]
    refpanel = file("https://github.com/ashotmarg/modules/raw/master/ashTestData/ref.22Jul22.46e.vcf.gz", checkIfExists: true)
    []
    []
    []

    BEAGLE_BEAGLE5 ( input, refpanel, [], [], [] )
}

workflow test_beagle_beagle5_ref_map {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/ashotmarg/modules/raw/master/ashTestData/target.22Jul22.46e.vcf.gz", checkIfExists: true)
    ]
    refpanel = file("https://github.com/ashotmarg/modules/raw/master/ashTestData/ref.22Jul22.46e.vcf.gz", checkIfExists: true)
    genmap   = file("https://github.com/ashotmarg/test-datasets/raw/modules/data/delete_me/beagle/plink.chr22.GRCh38.map", checkIfExists: true)
    []
    []

    BEAGLE_BEAGLE5 ( input, refpanel, genmap, [], [] )
}
