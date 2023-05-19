#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR          } from '../../../../../modules/nf-core/untar/main.nf'
include { CELLRANGERARC_MKFASTQ } from '../../../../../modules/nf-core/cellrangerarc/mkfastq/main.nf'

/*
Be ware! You do not have to split the cellranger-arc mkfastq call between ATAC and RNAseq files.
This totally depends on your experimental setup. Yet, the cellranger-arc datasets are split and
thus we have to call cellranger-arc mkfastq twice.
*/
workflow test_cellrangerarc_mkfastq_simple_atac {

    simple_csv_atac = file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-atac-simple-1.0.0.csv",
                            checkIfExists: true)
    tiny_bcl_atac = [ [], file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-atac-1.0.0.tar.gz",
                            checkIfExists: true) ]

    UNTAR ( tiny_bcl_atac )
    CELLRANGERARC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, simple_csv_atac)

}


workflow test_cellrangerarc_mkfastq_simple_gex {

    simple_csv_gex = file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-gex-simple-1.0.0.csv",
                            checkIfExists: true)
    tiny_bcl_gex = [ [], file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-gex-1.0.0.tar.gz",
                            checkIfExists: true) ]

    UNTAR ( tiny_bcl_gex )
    CELLRANGERARC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, simple_csv_gex)
}


workflow test_cellrangerarc_mkfastq_illumina_atac {

    samplesheet_csv_atac = file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-atac-samplesheet-1.0.0.csv",
                                checkIfExists: true)
    tiny_bcl_atac = [ [], file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-atac-1.0.0.tar.gz",
                                checkIfExists: true) ]

    UNTAR ( tiny_bcl_atac )
    CELLRANGERARC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, samplesheet_csv_atac)

}


workflow test_cellrangerarc_mkfastq_illumina_gex {

    samplesheet_csv_gex = file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-gex-samplesheet-1.0.0.csv",
                                checkIfExists: true)
    tiny_bcl_gex = [ [], file("https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-gex-1.0.0.tar.gz",
                                checkIfExists: true) ]

    UNTAR ( tiny_bcl_gex )
    CELLRANGERARC_MKFASTQ ( UNTAR.out.untar.map{ it[1] }, samplesheet_csv_gex)
}
