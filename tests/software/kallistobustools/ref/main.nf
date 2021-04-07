#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTOBUSTOOLS_REF } from '../../../../software/kallistobustools/ref/main.nf' addParams( options: [args:"standard"] )
include { KALLISTOBUSTOOLS_REF as KALLISTOBUSTOOLS_REF_LAMANNO } from '../../../../software/kallistobustools/ref/main.nf' addParams( options: [args:"lamanno"] )
include { KALLISTOBUSTOOLS_REF as KALLISTOBUSTOOLS_REF_NUCLEUS } from '../../../../software/kallistobustools/ref/main.nf' addParams( options: [args:"nucleus"] )

workflow test_kallistobustools_ref_standard {

    input = [ [id:'test.standard'], // meta map
        [file("${launchDir}/tests/data/delete_me/kallistobustools/GRCm39.chr19_100k.fa.gz", checkIfExists: true)]
        ]

    gtf      = file("${launchDir}/tests/data/delete_me/kallistobustools/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true) 
    
    KALLISTOBUSTOOLS_REF( input, gtf)
}

workflow test_kallistobustools_ref_lamanno {
    input = [ [id:'test.lamanno'], // meta map
              [file("${launchDir}/tests/data/delete_me/kallistobustools/GRCm39.chr19_100k.fa.gz", checkIfExists: true)]
            ]
    gtf   = file("${launchDir}/tests/data/delete_me/kallistobustools/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true) 
    KALLISTOBUSTOOLS_REF_LAMANNO( input, gtf)
}

workflow test_kallistobustools_ref_nucleus {

    input = [ [id:'test.nucleus'], // meta map
        [file("${launchDir}/tests/data/delete_me/kallistobustools/GRCm39.chr19_100k.fa.gz", checkIfExists: true)]
        ]

    gtf      = file("${launchDir}/tests/data/delete_me/kallistobustools/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true) 
    
    KALLISTOBUSTOOLS_REF_NUCLEUS( input, gtf)
}
