#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_EXTRACT  } from '../../../../../modules/nf-core/somalier/extract/main.nf'
include { SOMALIER_ANCESTRY } from '../../../../../modules/nf-core/somalier/ancestry/main.nf'
include { UNTAR             } from '../../../../../modules/nf-core/untar/main.nf'

workflow test_somalier_ancestry {

    labelled_somalier_tar = file("https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz")
    labels = file("https://github.com/brentp/somalier/raw/73db124d3fe9febe3a53787707554f863595b48f/scripts/ancestry-labels-1kg.tsv", checkIfExists: true)
    whole_genome_sommalier = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/whole_genome.somalier", checkIfExists: true)

    ch_labels = Channel.of([[id:"1kg"], labels])
    ch_whole_genome_somalier = Channel.of([[id:"test"], whole_genome_sommalier])

    UNTAR ( [[id:"1kg"], labelled_somalier_tar] )
    ch_labelled_somalier_files = ch_labels.join(UNTAR.out.untar)

    SOMALIER_ANCESTRY ( 
        ch_whole_genome_somalier, 
        ch_labelled_somalier_files 
    )

}
