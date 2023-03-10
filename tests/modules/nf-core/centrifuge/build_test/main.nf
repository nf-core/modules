#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CENTRIFUGE_BUILD } from '../../../../../modules/nf-core/centrifuge/build/main.nf'

workflow test_centrifuge_build {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file( 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
    ]

    conversion_table = file(params.test_data['sarscov2']['metagenome']['seqid2taxid_map'], checkIfExists: true)
    taxonomy_tree = file(params.test_data['sarscov2']['metagenome']['nodes_dmp'], checkIfExists: true)
    name_table = file(params.test_data['sarscov2']['metagenome']['names_dmp'], checkIfExists: true)

    CENTRIFUGE_BUILD ( input, conversion_table, taxonomy_tree, name_table )
}
