#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKM_LINEAGEWF              } from '../../../../modules/checkm/lineagewf/main.nf'
include { CHECKM_QA                     } from '../../../../modules/checkm/qa/main.nf'
include { CHECKM_QA as CHECKM_QA_FASTA  } from '../../../../modules/checkm/qa/main.nf'

workflow test_checkm_qa {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true) ]
    fasta_ext = 'fasta'

    CHECKM_LINEAGEWF ( input, fasta_ext, [] )

    ch_checkmqa_input = CHECKM_LINEAGEWF ( input, fasta_ext, [] )
        .join(CHECKM_LINEAGEWF.out.marker_file)
        .map{
            meta, dir, marker ->
            [ meta, dir, marker, []]
        }

    CHECKM_QA ( ch_checkmqa_input, [] )
}

workflow test_checkm_qa_fasta {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true) ]
    fasta_ext = 'fasta'

    ch_checkmqa_input = CHECKM_LINEAGEWF ( input, fasta_ext, [] )
        .join(CHECKM_LINEAGEWF.out.marker_file)
        .map{
            meta, dir, marker ->
            [ meta, dir, marker, []]
        }

    CHECKM_QA_FASTA ( ch_checkmqa_input, [] )
}
