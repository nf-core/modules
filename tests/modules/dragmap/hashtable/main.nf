#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DRAGMAP_HASHTABLE } from '../../../../modules/dragmap/hashtable/main.nf'

workflow test_dragmap_hashtable {
    
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    DRAGMAP_HASHTABLE ( fasta )
}

// TODO Add test using alt-masked bed file
// https://github.com/Illumina/dragmap#build-hash-table-using-an-alt-masked-bed-file
