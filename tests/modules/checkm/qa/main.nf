#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKM_LINEAGEWF              } from '../../../../modules/checkm/lineagewf/main.nf'
include { CHECKM_QA                     } from '../../../../modules/checkm/qa/main.nf'
include { CHECKM_QA as CHECKM_QA_FASTA  } from '../../../../modules/checkm/qa/main.nf'

workflow test_checkm_qa {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true) ]
    fasta_ext = 'fasta'

    CHECKM_LINEAGEWF ( input, fasta_ext, [] )
    CHECKM_QA ( CHECKM_LINEAGEWF.out.checkm_output, CHECK_LINEAGEWF.out.marker_file, [], [] )
}

workflow test_checkm_qa_fasta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true) ]
    fasta_ext = 'fasta'

    CHECKM_LINEAGEWF ( input, fasta_ext, [] )
    CHECKM_QA_FASTA ( CHECKM_LINEAGEWF.out.checkm_output, CHECK_LINEAGEWF.out.marker_file, [], [] )
}
