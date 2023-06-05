#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMOD_ANNOTATE } from '../../../../../modules/nf-core/genmod/annotate/main.nf'

input = [
    [ id:'test', single_end:false ], // meta map
    file(params.test_data['homo_sapiens']['illumina']['genmod_vcf_gz'], checkIfExists: true)
]

workflow test_genmod_annotate {

    GENMOD_ANNOTATE ( input )
}

workflow test_genmod_annotate_stub {

    GENMOD_ANNOTATE ( input )
}
