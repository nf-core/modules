#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMOD_SCORE    } from '../../../../../modules/nf-core/genmod/score/main.nf'

input = [
    [ id:'test', single_end:false ], // meta map
    file(params.test_data['homo_sapiens']['illumina']['genmod_models_vcf_gz'], checkIfExists: true)
]
fam = file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)
config = file(params.test_data['homo_sapiens']['illumina']['rank_model'], checkIfExists: true)

workflow test_genmod_score {

    GENMOD_SCORE    ( input, fam, config)
}

workflow test_genmod_score_stub {

    GENMOD_SCORE    ( input, fam, config)
}
