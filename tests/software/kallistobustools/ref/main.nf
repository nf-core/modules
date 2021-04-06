#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTOBUSTOOLS_REF } from '../../../../software/kallistobustools/ref/main.nf' addParams( options: [:] )

workflow test_kallistobustools_ref_standard {

    input = [ [id:'test.standard', workflow:"standard"], // meta map
        [file("${launchDir}/tests/data/delete_me/kallistobustools/GRCm39.chr19_100k.fa.gz", checkIfExists: true)]
        ]

    gtf      = file("${launchDir}/tests/data/delete_me/kallistobustools/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true) 
    
    KALLISTOBUSTOOLS_REF( input, gtf)
}

workflow test_kallistobustools_ref_lamanno {

    input = [ [id:'test.lamanno',workflow:"lamanno"], // meta map
        [file("${launchDir}/tests/data/delete_me/kallistobustools/GRCm39.chr19_100k.fa.gz", checkIfExists: true)]
        ]

    gtf      = file("${launchDir}/tests/data/delete_me/kallistobustools/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true) 
    
    KALLISTOBUSTOOLS_REF( input, gtf)
}

workflow test_kallistobustools_ref_nucleus {

    input = [ [id:'test.nucleus',workflow:"nucleus"], // meta map
        [file("${launchDir}/tests/data/delete_me/kallistobustools/GRCm39.chr19_100k.fa.gz", checkIfExists: true)]
        ]

    gtf      = file("${launchDir}/tests/data/delete_me/kallistobustools/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true) 
    
    KALLISTOBUSTOOLS_REF( input, gtf)
}
