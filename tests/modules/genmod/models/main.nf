#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMOD_ANNOTATE } from '../../../../modules/genmod/annotate/main.nf'
include { GENMOD_MODELS   } from '../../../../modules/genmod/models/main.nf'

workflow test_genmod_models {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['genmod_vcf_gz'], checkIfExists: true)
    ]
    fam = file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)

    GENMOD_ANNOTATE (input)
    GENMOD_MODELS ( GENMOD_ANNOTATE.out.vcf, fam, [] )
}

workflow test_genmod_models_stub {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['genmod_vcf_gz'], checkIfExists: true)
    ]
    fam = file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)

    GENMOD_ANNOTATE (input)
    GENMOD_MODELS ( GENMOD_ANNOTATE.out.vcf, fam, [] )
}
