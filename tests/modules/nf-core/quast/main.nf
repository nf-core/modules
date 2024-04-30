#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUAST as QUAST_REF         }   from '../../../../modules/nf-core/quast/main.nf'
include { QUAST as QUAST_NOREF_NOGFF }   from '../../../../modules/nf-core/quast/main.nf'
include { QUAST as QUAST_NOGFF       }   from '../../../../modules/nf-core/quast/main.nf'
include { QUAST as QUAST_NOREF       }   from '../../../../modules/nf-core/quast/main.nf'

workflow test_quast_ref {
    consensus = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)]
    fasta     = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)]
    gff       = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)]

    QUAST_REF ( consensus, fasta, gff )
}

workflow test_quast_noref_nogff {
    consensus = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)]
    fasta     = [[:],[]]
    gff       = [[:],[]]

    QUAST_NOREF_NOGFF ( consensus, fasta, gff )
}

workflow test_quast_nogff {
    consensus = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)]
    fasta     = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)]
    gff       = [[:],[]]

    QUAST_NOGFF ( consensus, fasta, gff )
}

workflow test_quast_noref {
    consensus = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)]
    fasta     = [[:],[]]
    gff       = [[ id:'test', single_end:false ], file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)]

    QUAST_NOREF ( consensus, fasta, gff )
}
