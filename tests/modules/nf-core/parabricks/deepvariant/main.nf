#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_DEEPVARIANT } from '../../../../../modules/nf-core/parabricks/deepvariant/main.nf'

workflow test_parabricks_deepvariant {
    
    input = [
        [ id:'test'],
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test.bam")
        
    ]

    PARABRICKS_DEEPVARIANT ( input, interval_file=[] )
}

workflow test_parabricks_deepvariant_intervals {
    
    input = [
        [ id:'test'],
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test.bam")
    ]
    interval_file = [
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/intervals_1.bed", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/intervals_2.bed", checkIfExists: true)
    ]

    PARABRICKS_DEEPVARIANT ( input, interval_file=interval_file )
}
