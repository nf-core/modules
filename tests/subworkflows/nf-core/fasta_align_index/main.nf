#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTA_ALIGN_INDEX } from '../../../../subworkflows/nf-core/fasta_align_index/main.nf'

workflow test_bowtie2_build {
    fasta = [
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    FASTA_ALIGN_INDEX ( fasta, 'bowtie2' )
}

workflow test_bwamem1_index {
    fasta = [
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    FASTA_ALIGN_INDEX ( fasta, 'bwamem' )
}

workflow test_bwamem2_index {
    fasta = [
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    FASTA_ALIGN_INDEX ( fasta, 'bwamem2' )
}

workflow test_dragmap_hashtable {
    fasta = [
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    FASTA_ALIGN_INDEX ( fasta, 'dragmap' )
}

workflow test_snap_index {
    fasta = [
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    FASTA_ALIGN_INDEX ( fasta, 'snap' )
}
