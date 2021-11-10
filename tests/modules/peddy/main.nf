#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PEDDY } from '../../../modules/peddy/main.nf' addParams( options: [:] )

workflow test_peddy {

    input = [ [ id:'test', single_end:false ],
            file('https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/data/genomics/homo_sapiens/genome/vcf/ped/justhusky.ped', checkIfExists: true)] // meta map
    vcf = file('https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/data/genomics/homo_sapiens/genome/vcf/ped/justhusky_minimal.vcf.gz', checkIfExists: true)
    tbi = file('https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/data/genomics/homo_sapiens/genome/vcf/ped/justhusky_minimal.vcf.gz.tbi', checkIfExists: true)

    PEDDY ( input , vcf , tbi )
}
