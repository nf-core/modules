#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SIMPLEAF_INDEX } from '../../../../../modules/nf-core/simpleaf/index/main.nf'

workflow test_simpleaf_index_expanded {

    genome_fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    meta = [ 'id': 'human'] 

    SIMPLEAF_INDEX ( 
        [ meta, genome_fasta ],
        [ meta, gtf ],
        [[],[]]
    )
}

workflow test_simpleaf_index_direct {

    transcriptome_fasta = file(params.test_data['homo_sapiens']['genome']['transcriptome_fasta'], checkIfExists: true)
    meta = [ 'id': 'human'] 

    SIMPLEAF_INDEX (
        [[],[]],
        [[],[]],
        [ meta, transcriptome_fasta ],
    )
}
