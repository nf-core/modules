#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUAST }   from '../../../modules/quast/main.nf'

workflow test_quast_ref {
    fasta     = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    gff       = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)
    consensus = file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)
    use_fasta = true
    use_gtf   = true

    QUAST ( consensus, fasta, gff, use_fasta, use_gtf )
}

workflow test_quast_noref {
    fasta     = file('fasta_dummy')
    gff       = file('gff_dummy')
    consensus = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    use_fasta = false
    use_gtf   = false

    QUAST ( consensus, fasta, gff, use_fasta, use_gtf )
}
