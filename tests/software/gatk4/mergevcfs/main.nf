#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MERGEVCFS } from '../../../../software/gatk4/mergevcfs/main.nf' addParams( options: [:] )

workflow test_gatk4_mergevcfs {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test2.vcf", checkIfExists: true) ] 
            ]
    dict = file('tests/data/genomics/sarscov2/genome/genome.dict', checkIfExists: true)

    GATK4_MERGEVCFS ( input, dict, false )
}

workflow test_gatk4_mergevcfs_refdict {
    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test2.vcf", checkIfExists: true) ] 
            ]
    dict  = file('tests/data/genomics/sarscov2/genome/genome.dict', checkIfExists: true)

    GATK4_MERGEVCFS ( input, ref_dict, true )
}
