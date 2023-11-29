#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMOD_MODELS   } from '../../../../../modules/nf-core/genmod/models/main.nf'

input = [
    [ id:'test', single_end:false ], // meta map
    file(params.test_data['homo_sapiens']['illumina']['genmod_annotate_vcf_gz'], checkIfExists: true)
]
fam = file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)

workflow test_genmod_models {

    GENMOD_MODELS ( input, fam, [] )
}

workflow test_genmod_models_stub {

    GENMOD_MODELS ( input, fam, [] )
}
