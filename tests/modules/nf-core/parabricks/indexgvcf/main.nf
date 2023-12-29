#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_INDEXGVCF } from '../../../../../modules/nf-core/parabricks/indexgvcf/main.nf'

workflow test_parabricks_indexgvcf {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
    ]

    process stage {
    input:
    tuple val(meta), path(gvcf)

    output:
    tuple val(meta), path("*.g.vcf")

    script:
    """
    mv $gvcf test.genome.g.vcf
    """
    }

    PARABRICKS_INDEXGVCF ( stage ( input ) )
}

workflow test_parabricks_indexgvcf_gz {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true)
    ]

    process stage_gz {
    stageInMode "copy"

    input:
    tuple val(meta), path(gvcf)

    output:
    tuple val(meta), path("*.g.vcf.gz")

    script:
    """
    mv $gvcf test.genome.g.vcf.gz
    """
    }

    PARABRICKS_INDEXGVCF ( stage_gz ( input ) )
}
