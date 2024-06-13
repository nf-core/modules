#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEAGLE5_BEAGLE } from '../../../../../modules/nf-core/beagle5/beagle/main.nf'

workflow test_beagle5_beagle {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/target.22Jul22.46e.vcf.gz", checkIfExists: true)
    ]
    refpanel    = []
    genmap      = []
    exclsamples = []
    exclmarkers = []

    BEAGLE5_BEAGLE ( input, refpanel, genmap, exclsamples, exclmarkers )
}

workflow test_beagle5_beagle_ref {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/target.22Jul22.46e.vcf.gz", checkIfExists: true)
    ]
    refpanel    = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/ref.22Jul22.46e.vcf.gz", checkIfExists: true)
    genmap      = []
    exclsamples = []
    exclmarkers = []

    BEAGLE5_BEAGLE ( input, refpanel, genmap, exclsamples, exclmarkers )
}

workflow test_beagle5_beagle_ref_map {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/target.22Jul22.46e.vcf.gz", checkIfExists: true)
    ]
    refpanel    = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/ref.22Jul22.46e.vcf.gz", checkIfExists: true)
    genmap      = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/beagle/plink.chr22.GRCh38.map", checkIfExists: true)
    exclsamples = []
    exclmarkers = []

    BEAGLE5_BEAGLE ( input, refpanel, genmap, exclsamples, exclmarkers )
}
