#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAPLOGREP2_CLASSIFY } from '../../../../../modules/nf-core/haplogrep2/classify/main.nf'

workflow test_haplogrep2_classify {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_mito_vcf'], checkIfExists: true)
    ]
    format = 'vcf'

    HAPLOGREP2_CLASSIFY ( input,format )
}
