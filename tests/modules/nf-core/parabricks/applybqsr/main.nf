#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_APPLYBQSR } from '../../../../../modules/nf-core/parabricks/applybqsr/main.nf'

workflow test_parabricks_applybqsr {
    
    input = [
        [ id:'test'],
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test.bam"),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test.table")
    ]

    PARABRICKS_APPLYBQSR ( input, interval_file=[] )
}

workflow test_parabricks_applybqsr_2 {
    
    input = [
        [ id:'test'],
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test.bam"),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.gz", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test.table")
    ]

    PARABRICKS_APPLYBQSR ( input, interval_file=[] )
}

workflow test_parabricks_applybqsr_intervals {
    
    input = [
        [ id:'test'],
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test.bam"),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test.table")
    ]
    interval_file = [
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/intervals_1.bed", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/intervals_2.bed", checkIfExists: true)
    ]
    PARABRICKS_APPLYBQSR ( input, interval_file=interval_file )
}
