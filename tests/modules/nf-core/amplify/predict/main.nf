#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRODIGAL } from '../../../../modules/nf-core/prodigal/main.nf'
include { AMPLIFY_PREDICT } from '../../../../../modules/nf-core/amplify/predict/main.nf'

workflow amplify_predict {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]
    model_dir = []

    PRODIGAL ( input, "gff" )
    AMPLIFY_PREDICT ( PRODIGAL.out.amino_acid_fasta, model_dir)
}
