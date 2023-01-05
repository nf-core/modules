#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEAGLE_BEAGLE5 } from '../../../../../modules/nf-core/beagle/beagle5/main.nf'

workflow test_beagle_beagle5 {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/target.22Jul22.46e.vcf.gz", checkIfExists: true)
    ]
    refpanel    = []
    genmap      = []
    exclsamples = []
    exclmarkers = []

    BEAGLE_BEAGLE5 ( input, refpanel, genmap, exclsamples, exclmarkers )
}

workflow test_beagle_beagle5_ref {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/target.22Jul22.46e.vcf.gz", checkIfExists: true)
    ]
    refpanel    = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/ref.22Jul22.46e.vcf.gz", checkIfExists: true)
    genmap      = []
    exclsamples = []
    exclmarkers = []

    BEAGLE_BEAGLE5 ( input, refpanel, genmap, exclsamples, exclmarkers )
}

workflow test_beagle_beagle5_ref_map {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/target.22Jul22.46e.vcf.gz", checkIfExists: true)
    ]
    refpanel    = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/ref.22Jul22.46e.vcf.gz", checkIfExists: true)
    genmap      = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/plink.chr22.GRCh38.map", checkIfExists: true)
    exclsamples = []
    exclmarkers = []

    BEAGLE_BEAGLE5 ( input, refpanel, genmap, exclsamples, exclmarkers )
}
