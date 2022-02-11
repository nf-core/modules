#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ICHORCNA_RUN } from '../../../../modules/ichorcna/run/main.nf'
include { HMMCOPY_READCOUNTER } from '../../../../modules/hmmcopy/readcounter/main.nf'
include { HMMCOPY_GCCOUNTER } from '../../../../modules/hmmcopy/gccounter/main.nf'
include { HMMCOPY_MAPCOUNTER } from '../../../../modules/hmmcopy/mapcounter/main.nf'
include { HMMCOPY_GENERATEMAP } from '../../../../modules/hmmcopy/generatemap/main.nf'

workflow test_ichorcna_run_no_panel {

    input = [ [ id:'test'], // meta map
              file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/MBC_315.ctDNA.reads.wig", checkIfExists: true)
            ]

    gcwig   = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/gc_hg19_1000kb.wig", checkIfExists: true)
    mapwig  = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/map_hg19_1000kb.wig", checkIfExists: true)

    panel_of_normals = []
    centromere = []

    ICHORCNA_RUN ( input, gcwig, mapwig, panel_of_normals, centromere)
}

workflow test_ichorcna_run_inc_panel {

    input = [ [ id:'test'], // meta map
              file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/MBC_315.ctDNA.reads.wig", checkIfExists: true)
            ]

    gcwig   = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/gc_hg19_1000kb.wig", checkIfExists: true)
    mapwig  = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/map_hg19_1000kb.wig", checkIfExists: true)

    panel_of_normals = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds", checkIfExists: true)

    centromere = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt", checkIfExists: true)

    ICHORCNA_RUN ( input, gcwig, mapwig, panel_of_normals, centromere)
}
