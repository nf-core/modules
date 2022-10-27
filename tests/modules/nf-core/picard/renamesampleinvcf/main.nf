#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { PICARD_RENAMESAMPLEINVCF } from "$moduleDir/modules/nf-core/picard/renamesampleinvcf/main.nf"

workflow test_picard_renamesampleinvcf {

    input = [ [ id:'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
            ]

    PICARD_RENAMESAMPLEINVCF ( input )
}
