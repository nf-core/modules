#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KB_REF } from '../../../../software/kb/ref/main.nf' addParams( options: [:] )

workflow test_kb_ref_standard {

    input = [ [id:'test.standard', workflow:"standard"], // meta map
        [file("${launchDir}/tests/data/delete_me/kb/GRCm39.chr19_100k.fa.gz", checkIfExists: true)]
        ]

    gtf      = file("${launchDir}/tests/data/delete_me/kb/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true) 
    
    KB_REF( input, gtf)
}

workflow test_kb_ref_lamanno {

    input = [ [id:'test.lamanno',workflow:"lamanno"], // meta map
        [file("${launchDir}/tests/data/delete_me/kb/GRCm39.chr19_100k.fa.gz", checkIfExists: true)]
        ]

    gtf      = file("${launchDir}/tests/data/delete_me/kb/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true) 
    
    KB_REF( input, gtf)
}

workflow test_kb_ref_nucleus {

    input = [ [id:'test.nucleus',workflow:"nucleus"], // meta map
        [file("${launchDir}/tests/data/delete_me/kb/GRCm39.chr19_100k.fa.gz", checkIfExists: true)]
        ]

    gtf      = file("${launchDir}/tests/data/delete_me/kb/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true) 
    
    KB_REF( input, gtf)
}
