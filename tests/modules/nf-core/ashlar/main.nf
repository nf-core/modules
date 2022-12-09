#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASHLAR } from '../../../../../modules/modules/nf-core/ashlar/main.nf'
// include { ASHLAR } from '../../../../../modules/modules/nf-core/ashlar/main.nf' addParams( options: [args: '--flip-mosaic-x'] )

def TEST_IMG = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-0{6,7}.ome.tiff"
def TEST_IMG_1 = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-06.ome.tiff"
def TEST_IMG_2 = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-07.ome.tiff"
def TEST_SHEET = '/home/pollen/github/modules/tests/modules/nf-core/ashlar/input_sheet.csv'

workflow test_ashlar {

    input_list =  [ [ [ id:'test', args: '--flip-y' ],
               "${TEST_IMG_1} ${TEST_IMG_2}" ] ]
    input_channel = Channel.fromList(input_list)

    ASHLAR ( input_channel )
}

include { INPUT_CHECK } from '../../../../../modules/modules/nf-core/ashlar/input_check.nf'

workflow test_ashlar_sheet {

    ch_input = file(TEST_SHEET)

    INPUT_CHECK (
        ch_input
    )
    .images
    .map {
        it -> [ [ id: it.sample, args: '--flip-y' ],
                it.file_list ]
    }
    .set { input_maps }

    ASHLAR ( input_maps )

}
