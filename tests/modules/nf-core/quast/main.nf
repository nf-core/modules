#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUAST }   from '../../../../modules/nf-core/quast/main.nf'

workflow test_quast_ref {
    consensus = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)]
    fasta     = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)]
    gff       = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)]

    QUAST ( consensus, fasta, gff )
}

workflow test_quast_noref {
    consensus = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)]
    fasta     = [[:],[]]
    gff       = [[:],[]]

    QUAST ( consensus, fasta, gff )
}
