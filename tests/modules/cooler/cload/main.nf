#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_CLOAD }                       from '../../../../modules/cooler/cload/main.nf' addParams( options: [args:'pairix'] )
include { COOLER_CLOAD as COOLER_CLOAD_PAIRS } from '../../../../modules/cooler/cload/main.nf' addParams( options: [args:'pairs --chrom1 1 --pos1 2 --chrom2 4 --pos2 5 -N'] )
include { COOLER_CLOAD as COOLER_CLOAD_TABIX } from '../../../../modules/cooler/cload/main.nf' addParams( options: [args:'tabix'] )


workflow test_cooler_cload_pairix {

    input = [ [ id:'test_pairix', single_end:false ], // meta map
             file(params.test_data['generic']['cooler']['test_pairix_pair_gz'], checkIfExists: true),
             file(params.test_data['generic']['cooler']['test_pairix_pair_gz_px2'], checkIfExists: true)]

    sizes    = file(params.test_data['generic']['cooler']['hg19_chrom_sizes'], checkIfExists: true)
    bin_size = 2000000

    COOLER_CLOAD ( input, bin_size, sizes )

}

workflow test_cooler_cload_pairs {

    input = [ [ id:'test_pairs', single_end:false ], // meta map
             file(params.test_data['generic']['cooler']['test_pairs_pair'], checkIfExists: true),
             []]

    sizes    = file(params.test_data['generic']['cooler']['hg19_chrom_sizes'], checkIfExists: true)
    bin_size = 2000000

    COOLER_CLOAD_PAIRS ( input, bin_size, sizes )

}

workflow test_cooler_cload_tabix {

    input = [ [ id:'test_tabix', single_end:false ], // meta map
             file(params.test_data['generic']['cooler']['test_tabix_pair_gz'], checkIfExists: true),
             file(params.test_data['generic']['cooler']['test_tabix_pair_gz_tbi'], checkIfExists: true)]

    sizes    = file(params.test_data['generic']['cooler']['hg19_chrom_sizes'], checkIfExists: true)
    bin_size = 2000000

    COOLER_CLOAD_TABIX ( input, bin_size, sizes )

}
