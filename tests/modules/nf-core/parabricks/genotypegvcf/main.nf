#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_GENOTYPEGVCF } from '../../../../../modules/nf-core/parabricks/genotypegvcf/main.nf'

workflow test_parabricks_genotypegvcf {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
    ]
    fasta = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    process stage {
        input:
        tuple val(meta), path(input)

        output:
        tuple val(meta), path("*.g.vcf")

        script:
        """
        mv $input test.genome.g.vcf
        """
    }

    PARABRICKS_GENOTYPEGVCF ( stage( input ), fasta )
}
