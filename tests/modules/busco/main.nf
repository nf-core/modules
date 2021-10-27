#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO } from '../../../modules/busco/main.nf' addParams( options: [args: '--mode genome --auto-lineage'] )

workflow test_busco {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file('/lustre/scratch123/tol/teams/team308/users/ps22/modules/genome.fa', checkIfExists: true) ]

    BUSCO ( input,
            [] )
}
