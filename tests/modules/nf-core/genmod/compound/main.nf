#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMOD_COMPOUND } from '../../../../../modules/nf-core/genmod/compound/main.nf'

input = [
    [ id:'test', single_end:false ], // meta map
    file(params.test_data['homo_sapiens']['illumina']['genmod_score_vcf_gz'], checkIfExists: true)
]
fam = file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)
config = file(params.test_data['homo_sapiens']['illumina']['rank_model'], checkIfExists: true)

workflow test_genmod_compound {

    GENMOD_COMPOUND ( input )
}

workflow test_genmod_compound_stub {

    GENMOD_COMPOUND ( input )
}
