#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { PRODIGAL } from "$moduleDir/modules/nf-core/prodigal/main.nf"
include { AMPLIFY_PREDICT } from "$moduleDir/modules/nf-core/amplify/predict/main.nf"

workflow amplify_predict {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]
    model_dir = []

    PRODIGAL ( input, "gff" )
    AMPLIFY_PREDICT ( PRODIGAL.out.amino_acid_fasta, model_dir)
}
