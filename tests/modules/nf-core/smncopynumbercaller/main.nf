#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SMNCOPYNUMBERCALLER } from '../../../../modules/nf-core/smncopynumbercaller/main.nf' addParams( options: [args: ''] )

workflow test_smncopynumbercaller {

    input = Channel.of([
        [ id:'test'], // meta map
        file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esG2AlnRep1.bam", checkIfExists: true),
        file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esG2AlnRep1.bam.bai", checkIfExists: true)
    ])

    SMNCOPYNUMBERCALLER ( input )
}
