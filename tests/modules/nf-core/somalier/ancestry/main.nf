#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_EXTRACT } from '../../../../../modules/nf-core/somalier/extract/main.nf'
include { SOMALIER_ANCESTRY } from '../../../../../modules/nf-core/somalier/ancestry/main.nf'
include { UNTAR      } from '../../../../modules/nf-core/untar/main.nf'

workflow test_somalier_ancestry {

    // Generate somalier files from bam
    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://storage.googleapis.com/genomics-public-data/simons-genome-diversity-project/vcf/LP6005441-DNA_A01.annotated.nh2.variants.vcf.gz", checkIfExists: true),
        file("https://storage.googleapis.com/genomics-public-data/simons-genome-diversity-project/vcf/LP6005441-DNA_A01.annotated.nh2.variants.vcf.gz.tbi", checkIfExists: true)
    ]

    fasta       = file("https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta", checkIfExists: true)
    fasta_fai   = file("https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai", checkIfExists: true)

    sites       = file("https://github.com/brentp/somalier/files/3412453/sites.hg19.vcf.gz", checkIfExists: true)

    SOMALIER_EXTRACT ( input, fasta, fasta_fai, sites )

    ch_query_somalier_files = SOMALIER_EXTRACT.out.extract

//    // Import reference labels and somalier files

    ch_labels_tsv = file("https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv", checkIfExists: true)

    labelled_somalier_files = [
        [],
        file("https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz", checkIfExists: true)
    ]

    UNTAR ( labelled_somalier_files )

    ch_labelled_somalier_files = UNTAR.out.untar


    UNTAR.out.untar.multiMap { meta, directory ->
        ch_meta_dir: [meta, directory]
        ch_directory: [directory]
    }
    .set { ch_untar_multimap}

    ch_labelled_somalier_files = ch_untar_multimap.ch_directory

    // Run somalier ancestry


    SOMALIER_ANCESTRY ( ch_query_somalier_files, ch_labels_tsv, ch_labelled_somalier_files )
}
