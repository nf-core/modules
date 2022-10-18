#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR          } from '../../../../../modules/nf-core/untar/main.nf'
include { UNIVERSC_MKFASTQ } from '../../../../../modules/nf-core/universc/cellranger_os_mkfastq/main.nf'

workflow test_universc_mkfastq_simple {

    simple_csv = file("https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv", checkIfExists: true)
    tiny_bcl = [ [], file("https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz", checkIfExists: true) ]

    UNTAR ( tiny_bcl )

    UNIVERSC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, simple_csv)
}

workflow test_universc_mkfastq_illumina {

    samplesheet_csv = file("https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-samplesheet-1.2.0.csv", checkIfExists: true)
    tiny_bcl = [ [], file("https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz", checkIfExists: true) ]

    UNTAR ( tiny_bcl )

    UNIVERSC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, samplesheet_csv)
}
