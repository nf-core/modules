#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_EXTRACT  } from '../../../../../modules/nf-core/somalier/extract/main.nf'
include { SOMALIER_ANCESTRY } from '../../../../../modules/nf-core/somalier/ancestry/main.nf'
include { UNTAR             } from '../../../../../modules/nf-core/untar/main.nf'

workflow test_somalier_ancestry {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://storage.googleapis.com/genomics-public-data/simons-genome-diversity-project/vcf/LP6005441-DNA_A01.annotated.nh2.variants.vcf.gz", checkIfExists: true),
        file("https://storage.googleapis.com/genomics-public-data/simons-genome-diversity-project/vcf/LP6005441-DNA_A01.annotated.nh2.variants.vcf.gz.tbi", checkIfExists: true)
    ]

    fasta       = file("https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta", checkIfExists: true)
    fasta_fai   = file("https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai", checkIfExists: true)

    sites       = file("https://github.com/brentp/somalier/files/3412453/sites.hg19.vcf.gz", checkIfExists: true)

    SOMALIER_EXTRACT ( input, fasta, fasta_fai, sites )

    // Import reference labels and somalier files
    labelled_somalier_tar = [
        [],
        file("https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz", checkIfExists: true)
    ]
    UNTAR ( labelled_somalier_tar )

    ch_labelled_somalier_files = [
        [id:'test'],
        file("https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv", checkIfExists: true),
        UNTAR.out.untar.map { meta, dir ->
            list = dir.eachFile {it}
            return list
        }
    ]

    // Run somalier ancestry
    SOMALIER_ANCESTRY ( SOMALIER_EXTRACT.out.extract, ch_labelled_somalier_files )
}
