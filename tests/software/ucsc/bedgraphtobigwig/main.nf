#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_BEDGRAPHTOBIGWIG  } from '../../../../software/ucsc/bedgraphtobigwig/main.nf' addParams( options: [:] )

workflow test_ucsc_bedgraphtobigwig {
    input = [ [ id:'test' ], // meta map
              [ file('https://raw.githubusercontent.com/igvteam/igv.js/master/test/data/wig/bedgraph-example-uscs.bedgraph', checkIfExists: true) ] 
            ]
    sizes = file('https://raw.githubusercontent.com/igvteam/igv.js/master/test/data/wig/chrom.sizes', checkIfExists: true)

    UCSC_BEDGRAPHTOBIGWIG ( input, sizes )
}
