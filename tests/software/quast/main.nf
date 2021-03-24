#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUAST }   from '../../../software/quast/main.nf'   addParams(options: [:])

workflow test_quast_ref {
    fasta     = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    gff       = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.gtf", checkIfExists: true)
    consensus = file("${launchDir}/tests/data/genomics/sarscov2/genome/transcriptome.fasta", checkIfExists: true)
    use_fasta = true
    use_gtf   = true

    QUAST ( consensus, fasta, gff, use_fasta, use_gtf )
}

workflow test_quast_noref {
    fasta     = file('fasta_dummy')
    gff       = file('gff_dummy')
    consensus = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    use_fasta = false
    use_gtf   = false

    QUAST ( consensus, fasta, gff, use_fasta, use_gtf )
}
