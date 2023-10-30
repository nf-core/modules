#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKM_LINEAGEWF } from '../../../../../modules/nf-core/checkm/lineagewf/main.nf'
include { CHECKM_QA        } from '../../../../../modules/nf-core/checkm/qa/main.nf'
include { GUNC_DOWNLOADDB  } from '../../../../../modules/nf-core/gunc/downloaddb/main.nf'
include { GUNC_MERGECHECKM } from '../../../../../modules/nf-core/gunc/mergecheckm/main.nf'
include { GUNC_RUN         } from '../../../../../modules/nf-core/gunc/run/main.nf'
workflow test_gunc_mergecheckm {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true) ]
    fasta_ext = 'fasta'

    // Run CheckM
    CHECKM_LINEAGEWF ( input, fasta_ext, [] )

    ch_checkmqa_input = CHECKM_LINEAGEWF.out.checkm_output
        .join(CHECKM_LINEAGEWF.out.marker_file)
        .map{
            meta, dir, marker ->
            [ meta, dir, marker, []]
        }

    CHECKM_QA ( ch_checkmqa_input, [] )

    // Run GUNC
    GUNC_DOWNLOADDB ( 'progenomes' )
    GUNC_RUN ( input, GUNC_DOWNLOADDB.out.db )

    // Merge

    ch_input_for_merge = GUNC_RUN.out.maxcss_level_tsv.join(CHECKM_QA.out.output)

    GUNC_MERGECHECKM ( ch_input_for_merge )
}
