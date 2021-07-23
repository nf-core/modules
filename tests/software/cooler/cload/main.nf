#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_CLOAD } from '../../../../software/cooler/cload/main.nf' addParams( options: [:] )

workflow test_cooler_cload {

    input = [ [ id:'test', single_end:false ], // meta map
             file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz", checkIfExists: true),
             file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz.px2", checkIfExists: true)]
    sizes = file(params.test_data['homo_sapiens']['genome']['genome_sizes'], checkIfExists: true)
    bin_size = 2000000
    COOLER_CLOAD ( input, bin_size, sizes )
}
