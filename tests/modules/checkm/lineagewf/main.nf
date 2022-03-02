#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKM_LINEAGEWF } from '../../../../modules/checkm/lineagewf/main.nf'

workflow test_checkm_lineagewf {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true) ]
    fasta_ext = 'fasta'

    CHECKM_LINEAGEWF ( input, fasta_ext )
}

workflow test_checkm_lineagewf_multi {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true),
                file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)] ]
    fasta_ext = 'fasta'

    CHECKM_LINEAGEWF ( input, fasta_ext )
}
