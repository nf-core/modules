#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFTOOLS as VCFTOOLS_BASE     } from '../../../software/vcftools/main.nf' addParams( options: ['args': '--freq']               )
include { VCFTOOLS as VCFTOOLS_OPTIONAL } from '../../../software/vcftools/main.nf' addParams( options: ['args': '--freq --exclude-bed'] )

workflow test_vcftools_vcf_base {
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf", checkIfExists: true) 
            ]

    VCFTOOLS_BASE ( input, [], [] )
}

workflow test_vcftools_vcfgz_base {
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf.gz", checkIfExists: true) 
            ]

    VCFTOOLS_BASE ( input, [], [] )
}

workflow test_vcftools_vcf_optional {
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf", checkIfExists: true) 
            ]
    bed   = file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true)

    VCFTOOLS_OPTIONAL ( input, bed, [] )
}

workflow test_vcftools_vcfgz_optional {
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf.gz", checkIfExists: true) 
            ]
    bed   = file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true)

    VCFTOOLS_OPTIONAL ( input, bed, [] )
}
