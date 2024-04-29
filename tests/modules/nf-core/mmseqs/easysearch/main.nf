#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MMSEQS_EASYSEARCH } from '../../../../../modules/nf-core/mmseqs/easysearch/main.nf'
include { MMSEQS_CREATEDB as MMSEQS_CREATEDB_TARGET } from '../../../../../modules/nf-core/mmseqs/createdb/main.nf'


workflow test_mmseqs_easysearch {
    
    input_query = [
        [ id:'test_query', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]
    

    input_target = [
        [ id:'test_target', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['scaffolds_fasta'], checkIfExists: true)
    ]
    ch_target_db = MMSEQS_CREATEDB_TARGET( input_target  ).db

    MMSEQS_EASYSEARCH ( input_query, ch_target_db )
}
