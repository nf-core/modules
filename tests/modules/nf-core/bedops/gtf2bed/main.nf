#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDOPS_GTF2BED } from '../../../../../modules/nf-core/bedops/gtf2bed/main.nf'

workflow test_bedops_gtf2bed {
    
    input = [
       //[ id:'test' ], // meta map
        file('https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/genes.gtf', checkIfExists: true)]

    BEDOPS_GTF2BED ( input )
}
