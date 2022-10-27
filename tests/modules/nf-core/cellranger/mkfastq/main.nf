#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { UNTAR          } from "$moduleDir/modules/nf-core/untar/main.nf"
include { CELLRANGER_MKFASTQ } from "$moduleDir/modules/nf-core/cellranger/mkfastq/main.nf"

workflow test_cellranger_mkfastq_simple {

    simple_csv = file("https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv", checkIfExists: true)
    tiny_bcl = [ [], file("https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz", checkIfExists: true) ]

    UNTAR ( tiny_bcl )

    CELLRANGER_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, simple_csv)
}

workflow test_cellranger_mkfastq_illumina {

    samplesheet_csv = file("https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-samplesheet-1.2.0.csv", checkIfExists: true)
    tiny_bcl = [ [], file("https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz", checkIfExists: true) ]

    UNTAR ( tiny_bcl )

    CELLRANGER_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, samplesheet_csv)
}
