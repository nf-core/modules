#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_SORTBAM } from '../../../../modules/fgbio/sortbam/main.nf'
include { FGBIO_CALLMOLECULARCONSENSUSREADS } from '../../../../modules/fgbio/callmolecularconsensusreads/main.nf'

workflow test_fgbio_callmolecularconsensusreads {
    input = [ [ id:'test' ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) ]
    FGBIO_SORTBAM ( input )
    FGBIO_CALLMOLECULARCONSENSUSREADS ( FGBIO_SORTBAM.out.bam  )
}
