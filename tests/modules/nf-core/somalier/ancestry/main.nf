#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_EXTRACT  } from '../../../../../modules/nf-core/somalier/extract/main.nf'
include { SOMALIER_ANCESTRY } from '../../../../../modules/nf-core/somalier/ancestry/main.nf'
include { UNTAR             } from '../../../../../modules/nf-core/untar/main.nf'

workflow test_somalier_ancestry {

    input = [
        [ id:'test', single_end:false ],
        file("https://storage.googleapis.com/genomics-public-data/simons-genome-diversity-project/vcf/LP6005441-DNA_A01.annotated.nh2.variants.vcf.gz", checkIfExists: true),
        file("https://storage.googleapis.com/genomics-public-data/simons-genome-diversity-project/vcf/LP6005441-DNA_A01.annotated.nh2.variants.vcf.gz.tbi", checkIfExists: true)
    ]

    fasta       = file("https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta", checkIfExists: true)
    fasta_fai   = file("https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai", checkIfExists: true)

    sites       = file("https://github.com/brentp/somalier/files/3412453/sites.hg19.vcf.gz", checkIfExists: true)

    labels      = file("https://github.com/brentp/somalier/raw/73db124d3fe9febe3a53787707554f863595b48f/scripts/ancestry-labels-1kg.tsv", checkIfExists: true)

    labelled_somalier_tar = [
        [id:"1kg"],
        file("https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz", checkIfExists: true)
    ]

    UNTAR ( labelled_somalier_tar )
    ch_labelled_somalier_files = labels.join(UNTAR.out.untar)

    SOMALIER_EXTRACT ( input, fasta, fasta_fai, sites )
    SOMALIER_ANCESTRY ( SOMALIER_EXTRACT.out.extract, ch_labelled_somalier_files )
}
