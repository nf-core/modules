#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAMTOOLS_SPLIT } from '../../../../modules/bamtools/split/main.nf' addParams( options: [args:"-reference"] )

workflow test_bamtools_split {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['alz1000'], checkIfExists: true) ]

    BAMTOOLS_SPLIT ( input )
}
