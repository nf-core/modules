#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR          } from '../../../../../modules/nf-core/untar/main.nf'
include { CELLRANGERARC_MKFASTQ } from '../../../../../modules/nf-core/cellrangerarc/mkfastq/main.nf'

workflow test_cellranger_arc_mkfastq_simple {

    simple_csv_atac = file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-atac-simple-1.0.0.csv",
                            checkIfExists: true)
    simple_csv_gex = file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-gex-simple-1.0.0.csv",
                            checkIfExists: true)
    tiny_bcl_atac = [ [], file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-atac-1.0.0.tar.gz",
                            checkIfExists: true) ]
    tiny_bcl_gex = [ [], file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-gex-1.0.0.tar.gz",
                            checkIfExists: true) ]

    UNTAR ( tiny_bcl_atac )
    CELLRANGERARC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, simple_csv_atac)

    UNTAR ( tiny_bcl_gex )
    CELLRANGERARC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, simple_csv_gex)
}

workflow test_cellranger_arc_mkfastq_illumina {

    samplesheet_csv_atac = file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-atac-samplesheet-1.0.0.csv",
                                checkIfExists: true)
    samplesheet_csv_gex = file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-gex-samplesheet-1.0.0.csv",
                                checkIfExists: true)
    tiny_bcl_atac = [ [], file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-atac-1.0.0.tar.gz",
                                checkIfExists: true) ]
    tiny_bcl_gex = [ [], file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-gex-1.0.0.tar.gz",
                                checkIfExists: true) ]

    UNTAR ( tiny_bcl_atac )
    CELLRANGERARC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, samplesheet_csv_atac)

    UNTAR ( tiny_bcl_gex )
    CELLRANGERARC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, samplesheet_csv_gex)
}
