#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTA_INDEX_DNA } from '../../../../subworkflows/nf-core/fasta_index_dna/main.nf'

workflow test_bowtie2_build {
    fasta = Channel.value([
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
    ])
    altliftover = Channel.value([
        [id:'test'],
        []
    ])
    FASTA_INDEX_DNA ( fasta, altliftover, 'bowtie2' )
}

workflow test_bwamem1_index {
    fasta = Channel.value([
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
    ])
    altliftover = Channel.value([
        [id:'test'],
        []
    ])
    FASTA_INDEX_DNA ( fasta, altliftover, 'bwamem' )
}

workflow test_bwamem2_index {
    fasta = Channel.value([
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
    ])
    altliftover = Channel.value([
        [id:'test'],
        []
    ])
    FASTA_INDEX_DNA ( fasta, altliftover, 'bwamem2' )
}

workflow test_dragmap_hashtable {
    fasta = Channel.value([
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
    ])
    altliftover = Channel.value([
        [id:'test'],
        []
    ])
    FASTA_INDEX_DNA ( fasta, altliftover, 'dragmap' )
}

workflow test_snap_index {
    fasta = Channel.value([
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
    ])
    altliftover = Channel.value([
        [id:'test'],
        []
    ])
    FASTA_INDEX_DNA ( fasta, altliftover, 'snap' )
}
