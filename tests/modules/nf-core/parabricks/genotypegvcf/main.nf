#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_GENOTYPEGVCF } from '../../../../../modules/nf-core/parabricks/genotypegvcf/main.nf'

workflow test_parabricks_genotypegvcf {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test.genome.vcf'], checkIfExists: true)
    ]
    fasta = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    PARABRICKS_GENOTYPEGVCF ( input, fasta )
}
