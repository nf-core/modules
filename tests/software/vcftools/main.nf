#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFTOOLS as VCFTOOLS_BASE} from '../../../software/vcftools/main.nf' addParams( options: ['args': '--freq'] )
include { VCFTOOLS as VCFTOOLS_OPTIONAL} from '../../../software/vcftools/main.nf' addParams( options: ['args': '--freq --exclude-bed'] )

workflow test_vcftools_vcf_base {

    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf", checkIfExists: true) ]

    VCFTOOLS_BASE ( input, [], [] )
}

workflow test_vcftools_vcfgz_base {

    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf.gz", checkIfExists: true) ]

    VCFTOOLS_BASE ( input, [], [] )
}

workflow test_vcftools_vcf_optional {

    def input = []
    def bed = file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true)

    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf", checkIfExists: true) ]

    VCFTOOLS_OPTIONAL ( input, bed, [] )
}

workflow test_vcftools_vcfgz_optional {

    def input = []
    def bed = file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true)

    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf.gz", checkIfExists: true) ]

    VCFTOOLS_OPTIONAL ( input, bed, [] )
}
