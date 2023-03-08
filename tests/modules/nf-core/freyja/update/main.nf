#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREYJA_UPDATE } from '../../../../../modules/nf-core/freyja/update/main.nf'

workflow test_freyja_update {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    FREYJA_UPDATE ( input )
}
