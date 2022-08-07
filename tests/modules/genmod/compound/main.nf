#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMOD_ANNOTATE } from '../../../../modules/genmod/annotate/main.nf'
include { GENMOD_MODELS   } from '../../../../modules/genmod/models/main.nf'
include { GENMOD_SCORE    } from '../../../../modules/genmod/score/main.nf'
include { GENMOD_COMPOUND } from '../../../../modules/genmod/compound/main.nf'

workflow test_genmod_compound {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['genmod_vcf_gz'], checkIfExists: true)
    ]
    fam = file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)
    config = file(params.test_data['homo_sapiens']['illumina']['rank_model'], checkIfExists: true)

    GENMOD_ANNOTATE (input)
    GENMOD_MODELS   ( GENMOD_ANNOTATE.out.vcf, fam, [] )
    GENMOD_SCORE    ( GENMOD_MODELS.out.vcf, fam, [], config)
    GENMOD_COMPOUND ( GENMOD_SCORE.out.vcf )
}

workflow test_genmod_compound_stub {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['genmod_vcf_gz'], checkIfExists: true)
    ]
    fam = file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)
    config = file(params.test_data['homo_sapiens']['illumina']['rank_model'], checkIfExists: true)

    GENMOD_ANNOTATE (input)
    GENMOD_MODELS   ( GENMOD_ANNOTATE.out.vcf, fam, [] )
    GENMOD_SCORE    ( GENMOD_MODELS.out.vcf, fam, [], config)
    GENMOD_COMPOUND ( GENMOD_SCORE.out.vcf )
}
