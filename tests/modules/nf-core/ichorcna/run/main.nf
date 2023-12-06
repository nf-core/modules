#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ICHORCNA_RUN } from '../../../../../modules/nf-core/ichorcna/run/main.nf'
include { HMMCOPY_READCOUNTER } from '../../../../../modules/nf-core/hmmcopy/readcounter/main.nf'
include { HMMCOPY_GCCOUNTER } from '../../../../../modules/nf-core/hmmcopy/gccounter/main.nf'
include { HMMCOPY_MAPCOUNTER } from '../../../../../modules/nf-core/hmmcopy/mapcounter/main.nf'
include { HMMCOPY_GENERATEMAP } from '../../../../../modules/nf-core/hmmcopy/generatemap/main.nf'

workflow test_ichorcna_run_no_panel {

    input = [ [ id:'test'], // meta map
              file("https://raw.githubusercontent.com/gavinhalab/ichorCNA/master/inst/extdata/MBC_315.ctDNA.reads.wig", checkIfExists: true)
            ]

    gc_wig   = file("https://raw.githubusercontent.com/gavinhalab/ichorCNA/master/inst/extdata/gc_hg19_1000kb.wig", checkIfExists: true)
    map_wig  = file("https://raw.githubusercontent.com/gavinhalab/ichorCNA/master/inst/extdata/map_hg19_1000kb.wig", checkIfExists: true)

    normal_wig       = []
    panel_of_normals = []
    rep_time_wig     = []
    exons            = []
    centromere       = file("https://raw.githubusercontent.com/gavinhalab/ichorCNA/master/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt", checkIfExists: true)

    ICHORCNA_RUN ( input, gc_wig, map_wig, normal_wig, panel_of_normals, rep_time_wig, exons, centromere )
}

workflow test_ichorcna_run_inc_panel {

    input = [ [ id:'test'], // meta map
              file("https://raw.githubusercontent.com/gavinhalab/ichorCNA/master/inst/extdata/MBC_315.ctDNA.reads.wig", checkIfExists: true)
            ]

    gc_wig   = file("https://raw.githubusercontent.com/gavinhalab/ichorCNA/master/inst/extdata/gc_hg19_1000kb.wig", checkIfExists: true)
    map_wig  = file("https://raw.githubusercontent.com/gavinhalab/ichorCNA/master/inst/extdata/map_hg19_1000kb.wig", checkIfExists: true)

    normal_wig       = []
    panel_of_normals = file("https://raw.githubusercontent.com/gavinhalab/ichorCNA/master/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds", checkIfExists: true)
    rep_time_wig     = []
    exons            = []
    centromere = file("https://raw.githubusercontent.com/gavinhalab/ichorCNA/master/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt", checkIfExists: true)

    ICHORCNA_RUN ( input, gc_wig, map_wig, normal_wig, panel_of_normals, rep_time_wig, exons, centromere)
}
