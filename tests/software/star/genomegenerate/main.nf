#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../../../../software/star/genomegenerate/main.nf' addParams( options: [args: '--genomeSAindexNbases 9'] )

workflow test_star_genomegenerate {
    fasta = file("${launchDir}/tests/data/generic/fasta/GCF_000019425.1_ASM1942v1_genomic.fna", checkIfExists: true)
    gtf   = file("${launchDir}/tests/data/generic/gtf/GCF_000019425.1_ASM1942v1_genomic.gtf", checkIfExists: true)

    STAR_GENOMEGENERATE ( fasta, gtf )
}
