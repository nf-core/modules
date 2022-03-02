#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QCAT } from '../../../modules/qcat/main.nf'

workflow test_qcat {
    input = [ [ id:'test' ], // meta map
              [ file("https://github.com/nf-core/test-datasets/raw/nanoseq/fastq/nondemultiplexed/sample_nobc_dx.fastq.gz", checkIfExists: true) ]
            ]
    barcode_kit = 'NBD103/NBD104'

    QCAT ( input, barcode_kit )
}
