#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ICHORCNA_CREATEPON } from '../../../../modules/ichorcna/createpon/main.nf'

workflow test_ichorcna_createpon {

    input = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/MBC_315.ctDNA.reads.wig", checkIfExists: true)

    gcwig   = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/gc_hg19_1000kb.wig", checkIfExists: true)
    mapwig  = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/map_hg19_1000kb.wig", checkIfExists: true)

    centromere = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt", checkIfExists: true)

    ICHORCNA_CREATEPON ( input, gcwig, mapwig, centromere )
}

workflow test_ichorcna_createpon2 {

    input = [file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/MBC_315.ctDNA.reads.wig", checkIfExists: true),
             file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/MBC_315_T2.ctDNA.reads.wig", checkIfExists: true)]

    gcwig   = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/gc_hg19_1000kb.wig", checkIfExists: true)
    mapwig  = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/map_hg19_1000kb.wig", checkIfExists: true)

    centromere = file("https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt", checkIfExists: true)

    ICHORCNA_CREATEPON ( input, gcwig, mapwig, centromere )
}
