#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../../../../software/star/genomegenerate/main.nf' addParams( options: [args: '--genomeSAindexNbases 9'] )

workflow test_star_genomegenerate {

    def fasta = file("${launchDir}/tests/data/fasta/E_coli/GCF_000019425.1_ASM1942v1_genomic.fna", checkIfExists: true)
    def gtf   = file("${launchDir}/tests/data/gff/GCF_000019425.1_ASM1942v1_genomic.gtf", checkIfExists: true)
    STAR_GENOMEGENERATE ( fasta, gtf )
}
