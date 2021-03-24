#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

test_options = ['args': '--filter-name "test_filter" --filter-expression "MQ0 > 0"']
include { GATK4_VARIANTFILTRATION } from '../../../../software/gatk4/variantfiltration/main.nf' addParams( options: test_options )

workflow test_gatk4_variantfiltration {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf", checkIfExists: true) ] 
            ]
    fasta = [ file("tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true),
              file("tests/data/genomics/sarscov2/genome/genome.fasta.fai", checkIfExists: true),
              file("tests/data/genomics/sarscov2/genome/genome.dict", checkIfExists: true)
            ]

    GATK4_VARIANTFILTRATION ( input, fasta )
}
