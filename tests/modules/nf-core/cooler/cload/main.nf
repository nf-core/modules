#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_CLOAD }                       from '../../../../../modules/nf-core/cooler/cload/main.nf'
include { COOLER_CLOAD as COOLER_CLOAD_PAIRS } from '../../../../../modules/nf-core/cooler/cload/main.nf'
include { COOLER_CLOAD as COOLER_CLOAD_TABIX } from '../../../../../modules/nf-core/cooler/cload/main.nf'
include { COOLER_DUMP }                        from '../../../../../modules/nf-core/cooler/dump/main.nf'
include { COOLER_DUMP as COOLER_DUMP_PAIRS}    from '../../../../../modules/nf-core/cooler/dump/main.nf'
include { COOLER_DUMP as COOLER_DUMP_TABIX}    from '../../../../../modules/nf-core/cooler/dump/main.nf'

workflow test_cooler_cload_pairix {

    input = [ [ id:'test_pairix', single_end:false ], // meta map
             file(params.test_data['generic']['cooler']['test_pairix_pair_gz'], checkIfExists: true),
             file(params.test_data['generic']['cooler']['test_pairix_pair_gz_px2'], checkIfExists: true),
             2000000
    ]

    sizes    = file(params.test_data['generic']['cooler']['hg19_chrom_sizes'], checkIfExists: true)

    COOLER_CLOAD ( input, sizes )
    COOLER_DUMP(COOLER_CLOAD.out.cool.map { [it[0], it[1], []] })

}

workflow test_cooler_cload_pairs {

    input = [ [ id:'test_pairs', single_end:false ], // meta map
             file(params.test_data['generic']['cooler']['test_pairs_pair'], checkIfExists: true),
             [],
             2000000
    ]

    sizes    = file(params.test_data['generic']['cooler']['hg19_chrom_sizes'], checkIfExists: true)

    COOLER_CLOAD_PAIRS ( input, sizes )
    COOLER_DUMP_PAIRS(COOLER_CLOAD_PAIRS.out.cool.map { [it[0], it[1], []] })

}

workflow test_cooler_cload_tabix {

    input = [ [ id:'test_tabix', single_end:false ], // meta map
             file(params.test_data['generic']['cooler']['test_tabix_pair_gz'], checkIfExists: true),
             file(params.test_data['generic']['cooler']['test_tabix_pair_gz_tbi'], checkIfExists: true),
             2000000
    ]

    sizes    = file(params.test_data['generic']['cooler']['hg19_chrom_sizes'], checkIfExists: true)

    COOLER_CLOAD_TABIX ( input, sizes )
    COOLER_DUMP_TABIX(COOLER_CLOAD_TABIX.out.cool.map { [it[0], it[1], []] })

}
