#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HUMANN } from '../../../modules/humann/main.nf'

workflow test_humann {
    
    input = [
    	      [ id:'mpa_v30_CHOCOPhlAn_201901' ], // meta map
              file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/umi-dna/qiaseq/SRR7545951-small_1.fastq.gz', checkIfExists: true),
              file('chocophlan', checkIfExists: true),
              file('uniref', checkIfExists: true),
              file('metaphlandb', checkIfExists: true)
    ]

    HUMANN ( input )
}
