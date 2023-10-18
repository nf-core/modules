#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_INDEXGVCF } from '../../../../../modules/nf-core/parabricks/indexgvcf/main.nf'

workflow test_parabricks_indexgvcf {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
    ]

    PARABRICKS_INDEXGVCF ( input )
}

workflow test_parabricks_indexgvcf_gz {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true)
    ]

    PARABRICKS_INDEXGVCF ( input )
}
