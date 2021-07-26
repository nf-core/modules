#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_CLOAD } from '../../../../modules/cooler/cload/main.nf' addParams( options: [:] )
include { COOLER_CLOAD
 as COOLER_CLOAD_HICLIB} from '../../../../modules/cooler/cload/main.nf' addParams( options: [args:'hiclib'] )
include { COOLER_CLOAD
 as COOLER_CLOAD_PAIRS } from '../../../../modules/cooler/cload/main.nf' addParams( options: [args:'pairs --chrom1 1 --pos1 1 --chrom2 1 --pos2 1'] )
include { COOLER_CLOAD
 as COOLER_CLOAD_TABIX } from '../../../../modules/cooler/cload/main.nf' addParams( options: [args:'tabix'] )


workflow test_cooler_cload {

    input = [ [ id:'test', single_end:false ], // meta map
             file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz", checkIfExists: true),
             file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz.px2", checkIfExists: true)]
    sizes = file(params.test_data['homo_sapiens']['genome']['genome_sizes'], checkIfExists: true)
    bin_size = 2000000
    COOLER_CLOAD ( input, bin_size, sizes )

    input2 = [ [ id:'test', single_end:false ], // meta map
             file("https://raw.githubusercontent.com/mirnylab/hiclib-legacy/master/tests/fragmentHiC/test-hg19.hdf5", checkIfExists: true),
             []]
    //COOLER_CLOAD_HICLIB ( input2, bin_size, sizes )

    input3 = [ [ id:'test', single_end:false ], // meta map
             file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/hg19.sample1.pairs", checkIfExists: true),
             []]
    //COOLER_CLOAD_PAIRS ( input3, bin_size, sizes )

    input4 = [ [ id:'test', single_end:false ], // meta map
             file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/hg19.GM12878-MboI.pairs.subsample.sorted.txt.gz", checkIfExists: true),
             file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/hg19.GM12878-MboI.pairs.subsample.sorted.txt.gz.tbi", checkIfExists: true)]
    COOLER_CLOAD_TABIX ( input4, bin_size, sizes )

}
