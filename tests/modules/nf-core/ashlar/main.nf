#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// optional ags go below (if needed): https://nf-co.re/docs/contributing/tutorials/dsl2_modules_tutorial#passing-options-args
// TO DO: if they apply to all runs of this workflow they can't they also be in ext.args in nextflow.config?
include { ASHLAR } from '../../../../../modules/modules/nf-core/ashlar/main.nf' addParams( options: [args: '--flip-mosaic-x'] )
// include { ASHLAR } from '../../../../../modules/modules/nf-core/ashlar/main.nf'

def TEST_IMG = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-0{6,7}.ome.tiff"
// def TEST_IMG = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-06.ome.tiff"
def TEST_SHEET = '/home/pollen/github/modules/tests/modules/nf-core/ashlar/input_sheet.csv'

workflow test_ashlar {

    input_list =  [ [ [ id:'test', args: '--flip-y' ],
               file(TEST_IMG, checkIfExists: true) ] ]
    input_channel = Channel.fromList(input_list)

    ASHLAR ( input_channel )
}


// we can add additional test workflows below

include { INPUT_CHECK } from '../../../../../modules/modules/nf-core/ashlar/input_check.nf'

workflow test_ashlar_sheet {

    ch_input = file(TEST_SHEET)

    INPUT_CHECK (
        ch_input
    )
    .images
    .map {
        [ [ id:it.id, args: it.args],
            file(it.file_list, checkIfExists: true) ]
    }
    .set { input_maps }

    ASHLAR ( input_maps )

}

