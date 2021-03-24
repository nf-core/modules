#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OPTITYPE } from ''../../../software/optitype/main.nf' addParams( options: ['args':'-e 1 -b 0.009', 'args2':'solver=glpk'] )

workflow test_optitype {
    input = [ [ id:'test', seq_type:'dna' ], // meta map
              file("${launchDir}/tests/data/genomics/homo_sapiens/illumina/bam/example_pe.bam", checkIfExists: true) 
            ]
            
    OPTITYPE ( input )
}
