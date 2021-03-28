#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_SORTBAM } from '../../../../software/fgbio/sortbam/main.nf' addParams( options: [args: '-s TemplateCoordinate', suffix: '_out'] )
include { FGBIO_CALLMOLECULARCONSENSUSREADS } from '../../../../software/fgbio/callmolecularconsensusreads/main.nf' addParams( options: [args: '-M 1', suffix: '_molreads'] )

workflow test_fgbio_callmolecularconsensusreads {
    input = [ [ id:'test' ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) ]
    FGBIO_SORTBAM ( input )
    FGBIO_CALLMOLECULARCONSENSUSREADS ( FGBIO_SORTBAM.out.bam  )
}
