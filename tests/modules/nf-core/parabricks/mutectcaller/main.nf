#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_MUTECTCALLER } from '../../../../../modules/nf-core/parabricks/mutectcaller/main.nf'

workflow test_parabricks_mutectcaller {
    
    input = [
        [ id:'test' ],
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test.bam")        
    ]
    input2 = [
        [],
        []
    ]
    PARABRICKS_MUTECTCALLER ( input, input2=input2 interval_file=[] )
}

workflow test_parabricks_mutectcaller_tn {
    
    input = [
        [ id:'test' ],
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test.bam")
    ]
    input2 = [
        [ id:'test2' ],
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test2.bam")
    ]

    PARABRICKS_MUTECTCALLER ( input, input2=input2, interval_file=[] )
}


workflow test_parabricks_mutectcaller_intervals {
    
    input = [
        [ id:'test' ],
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/nf-core-files/test.bam")
    ]
    input2 = [
        [],
        []
    ]
    interval_file = [
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/intervals_1.bed", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/intervals_2.bed", checkIfExists: true)
    ]

    PARABRICKS_MUTECTCALLER ( input, input2=input2, interval_file=interval_file )
}
