#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QCAT  } from '../../../software/qcat/main.nf' addParams( options: [:] )

workflow test_qcat {
    def input = []
    input = [ [ id:'test' ], // meta map
            [ file(params.test_data['homo_sapiens']['nanopore']['non_demultiplexed_fastq'], checkIfExists: true) ]
            ]
    barcode_kit = 'NBD103/NBD104'
    QCAT ( input, barcode_kit )
}
