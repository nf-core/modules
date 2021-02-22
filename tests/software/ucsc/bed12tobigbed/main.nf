#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_BED12TOBIGBED } from '../../../../software/ucsc/bed12tobigbed/main.nf' addParams( options: [:] )

workflow test_ucsc_bed12tobigbed {

    def input = []
    input = [ 'test' , file("${launchDir}/tests/data/chrom_size/hg19.chrom.sizes", checkIfExists: true),
              'FALSE' , file("${launchDir}/tests/data/bed/test_ucsc_bed12tobigbed.bed", checkIfExists: true)]

    UCSC_BED12TOBIGBED ( input )
}
