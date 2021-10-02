#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GTDBTK_CLASSIFYWF } from '../../../../modules/gtdbtk/classifywf/main.nf' addParams( options: [:] )

process STUB_GTDBTK {
    output:
    path("*fa")                                                   , emit: bins
    tuple val("gtdbtk_r202_data"), path("database/*")             , emit: database

    stub:
    """
    touch 1.fa 2.fa 3.fa

    mkdir database
    touch database/gtdbtk_r202_data
    """
}

workflow test_gtdbtk_classifywf {

    STUB_GTDBTK()
    
    input = [[ id:'test', single_end:false, assembler:'SPADES' ], // meta map
             STUB_GTDBTK.out.bins.groupTuple()]

    GTDBTK_CLASSIFYWF ( input, STUB_GTDBTK.out.database )
}
