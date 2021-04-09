#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTOBUSTOOLS_REF } from '../../../../software/kallistobustools/ref/main.nf' addParams( options: [:] )

workflow test_kallistobustools_ref_standard {

    fasta       = file("${launchDir}/tests/data/delete_me/kallistobustools/GRCm39.chr19_100k.fa.gz", checkIfExists: true)
    gtf         = file("${launchDir}/tests/data/delete_me/kallistobustools/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true)
    workflow    = "standard"
    
    KALLISTOBUSTOOLS_REF(fasta, gtf, workflow)
}

workflow test_kallistobustools_ref_lamanno {

    fasta       = file("${launchDir}/tests/data/delete_me/kallistobustools/GRCm39.chr19_100k.fa.gz", checkIfExists: true)
    gtf         = file("${launchDir}/tests/data/delete_me/kallistobustools/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true)
    workflow    = "standard"
    
    KALLISTOBUSTOOLS_REF( fasta, gtf, workflow)
}

workflow test_kallistobustools_ref_nucleus {

    fasta       = file("${launchDir}/tests/data/delete_me/kallistobustools/GRCm39.chr19_100k.fa.gz", checkIfExists: true)
    gtf         = file("${launchDir}/tests/data/delete_me/kallistobustools/gencode.VM26.chr19_10k.gtf.gz", checkIfExists: true)
    workflow    = "standard"
    
    KALLISTOBUSTOOLS_REF( fasta, gtf, workflow)
}
