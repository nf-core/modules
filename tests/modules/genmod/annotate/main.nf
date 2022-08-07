#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMOD_ANNOTATE } from '../../../../modules/genmod/annotate/main.nf'

workflow test_genmod_annotate {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['genmod_vcf_gz'], checkIfExists: true)
    ]

    GENMOD_ANNOTATE ( input )
}

workflow test_genmod_annotate_stub {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['genmod_vcf_gz'], checkIfExists: true)
    ]

    GENMOD_ANNOTATE ( input )
}
