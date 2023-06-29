#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR          } from '../../../../../modules/nf-core/untar/main.nf'
include { CELLRANGERATAC_MKFASTQ } from '../../../../../modules/nf-core/cellrangeratac/mkfastq/main.nf'

workflow test_cellrangeratac_mkfastq_simple {

    simple_csv = file("https://cf.10xgenomics.com/supp/cell-atac/cellranger-atac-tiny-bcl-simple-1.0.0.csv", checkIfExists: true)
    tiny_bcl = [ [], file("https://cf.10xgenomics.com/supp/cell-atac/cellranger-atac-tiny-bcl-1.0.0.tar.gz", checkIfExists: true) ]

    UNTAR ( tiny_bcl )

    CELLRANGERATAC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, simple_csv)
}

workflow test_cellrangeratac_mkfastq_illumina {

    samplesheet_csv = file("https://cf.10xgenomics.com/supp/cell-atac/cellranger-atac-tiny-bcl-samplesheet-1.0.0.csv", checkIfExists: true)
    tiny_bcl = [ [], file("https://cf.10xgenomics.com/supp/cell-atac/cellranger-atac-tiny-bcl-1.0.0.tar.gz", checkIfExists: true) ]

    UNTAR ( tiny_bcl )

    CELLRANGERATAC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, samplesheet_csv)
}
