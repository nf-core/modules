#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CIRCEXPLORER2_ANNOTATE } from '../../../../../modules/nf-core/circexplorer2/annotate/main.nf'

workflow test_circexplorer2_annotate {

    input = [
        [ id:'fust1_3' ],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/circexplorer2/fust1_3.junction.bed")
    ]

    fasta = [
        file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/reference/chrI.fa")
    ]

    gene_annotation = [
        file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/reference/chrI.txt")
    ]

    CIRCEXPLORER2_ANNOTATE ( input, fasta, gene_annotation )
}
